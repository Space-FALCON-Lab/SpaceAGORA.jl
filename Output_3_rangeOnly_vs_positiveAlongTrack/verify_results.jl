import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration")))

using Arrow, CSV, DataFrames, Dates, LinearAlgebra, Test, TOML
include("simulation.jl")
using .ComparisonSimulation
const CS = ComparisonSimulation
const ROOT = @__DIR__
const VALIDATION = joinpath(ROOT, "validation")
mkpath(VALIDATION)

refinement = NamedTuple[]
settings = Settings()
expected_initial = initial_state(settings)
fine_settings = Settings(max_step_s=0.5)

@testset "Saved 80-orbit comparison and full-duration refinement" begin
    for policy in (:naive_next_entering, :positive_along_track)
        @testset "$policy" begin
            directory = joinpath(ROOT, "results", string(policy))
            info = TOML.parsefile(joinpath(directory, "run_info.toml"))
            @test info["status"] == "verified"
            @test info["settings"] == settings_dict(settings)
            @test info["settings"]["helpers"] == 200
            @test info["settings"]["orbits"] == 80
            @test info["settings"]["j2"] == J2
            @test info["settings"]["target_altitude_km"] == 1000
            @test info["settings"]["helper_altitude_km"] == 1050
            @test info["settings"]["target_inclination_deg"] == 0.5
            @test info["settings"]["helper_inclination_deg"] == 0
            @test info["trajectory_rows"] == 16001
            @test info["trajectory_columns"] == 1208
            for (name, checksum) in info["checksums"]
                @test CS.filehash(joinpath(directory, name)) == checksum
            end
            for (name, checksum) in info["source_sha256"]
                @test CS.filehash(joinpath(ROOT, name)) == checksum
            end
            trajectory = Arrow.Table(joinpath(directory, "trajectory.feather"))
            table = CSV.read(joinpath(directory, "target_orbital_elements.csv"), DataFrame)
            @test table.time_s == trajectory.time_s
            @test table.target_orbits == trajectory.target_orbits
            @test first(table.time_s) == 0
            @test last(table.time_s) == duration(settings)
            @test last(table.target_orbits) ≈ 80 atol=1e-12 rtol=0
            for index in 1:201, (field, offset) in (("pos", 0), ("vel", 3)), axis in 1:3
                column = trajectory[Symbol("sc$(index)_$(field)_$axis")]
                @test first(column) == expected_initial[offset + axis, index]
                @test all(isfinite, column)
            end
            for row in 1:nrow(table)
                pos = [trajectory[Symbol("sc1_pos_$axis")][row] for axis in 1:3]
                vel = [trajectory[Symbol("sc1_vel_$axis")][row] for axis in 1:3]
                oe = orbital_elements(pos, vel)
                @test oe.a_km ≈ table.a_km[row] atol=1e-10 rtol=0
                @test oe.eccentricity ≈ table.eccentricity[row] atol=1e-14 rtol=0
                @test oe.inclination_deg ≈ table.inclination_deg[row] atol=1e-12 rtol=0
                angle_error = deg2rad(oe.raan_deg - table.raan_deg[row])
                @test abs(atan(sin(angle_error), cos(angle_error))) < 1e-11
            end
            history = CSV.read(joinpath(directory, "schedule.csv"), DataFrame)
            @test nrow(history) == info["accepted_steps"] + 1
            @test all(>(0), diff(history.time_s))
            selected = history.active_helper_idx .> 0
            @test all(x -> 2 <= x <= 201, history.active_helper_idx[selected])
            @test all(x -> 0 <= x <= settings.range_m, history.distance_m[selected])
            if policy === :positive_along_track
                @test all(>(0), history.along_track_projection[selected])
            end

            fine = propagate(fine_settings, policy)
            fine_table = target_table(fine)
            @test fine_table.time_s == table.time_s
            CSV.write(joinpath(VALIDATION, "$(policy)_0p5s_target_elements.csv"), fine_table)
            differences = (max_a_difference_m=1000maximum(abs, table.a_km - fine_table.a_km),
                max_e_difference=maximum(abs, table.eccentricity - fine_table.eccentricity),
                max_i_difference_deg=maximum(abs, table.inclination_deg - fine_table.inclination_deg),
                max_raan_difference_deg=maximum(abs, table.raan_deg - fine_table.raan_deg),
                final_a_difference_m=1000(table.a_km[end] - fine_table.a_km[end]))
            push!(refinement, merge((policy=string(policy), orbits=80.0,
                coarse_max_step_s=settings.max_step_s, fine_max_step_s=fine_settings.max_step_s,
                refined_accepted_steps=fine.solution.stats.naccept), differences))
            println("$policy full-duration refinement: $differences")
            @test differences.max_a_difference_m < 10.0
            @test differences.max_e_difference < 2e-6
            @test differences.max_i_difference_deg < 1e-4
            @test differences.max_raan_difference_deg < 0.01
        end
    end
end

CSV.write(joinpath(VALIDATION, "full_duration_step_refinement.csv"), DataFrame(refinement))
CS.write_info(joinpath(VALIDATION, "production_verification.toml"),
    Dict("passed" => true, "finished_utc" => string(now(UTC)),
        "policies" => collect(string.(CS.POLICIES)), "verified_rows_per_policy" => 16001,
        "verified_spacecraft_per_policy" => 201,
        "step_refinement" => "Full 80 orbits; 1 s versus 0.5 s maximum accepted step",
        "comparison_thresholds" => Dict("a_m" => 10.0, "eccentricity" => 2e-6,
            "inclination_deg" => 1e-4, "raan_deg" => 0.01),
        "scope_note" => "Empirical step-refinement comparison, not a rigorous absolute error bound",
        "verification_source_sha256" => CS.filehash(@__FILE__)))
