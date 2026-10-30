module PeriapsisCorridorRuntimeTests

using Test
using CSV
using DataFrames
using TOML

# Keep the example definitions in a module so the script's PROGRAM_FILE guard
# cannot launch a mission while this test file is included by runtests.jl.
include(joinpath(@__DIR__, "../../../examples/AGORA_Earth_Periapsis_Corridor.jl"))

function _fresh_corridor_config(; frequent=false)
    source = load_corridor_config()
    values = deepcopy(source.values)
    if frequent
        values["comparison"]["frequent_deadband_m"] = 0.1
        values["comparison"]["target_periapsis_altitude_m"] = 98_500.0
    end
    config_root = mktempdir()
    config_path = joinpath(config_root, "periapsis_corridor_config.toml")
    open(config_path, "w") do io
        TOML.print(io, values)
    end
    return load_corridor_config(config_path)
end

function _case_output(root::AbstractString, name::AbstractString)
    return joinpath(root, name)
end

@testset "physical corridor finite-burn runtime" begin
    @testset "frequent policy delivers separate finite burns with logged mass and Δv" begin
        config = _fresh_corridor_config(frequent=true)

        apoapsis_m = 6_378_136.6 + Float64(config.values["orbit"]["apoapsis_altitude_m"])
        periapsis_m = 6_378_136.6 + Float64(config.values["orbit"]["initial_periapsis_altitude_m"])
        semimajor_axis_m = (apoapsis_m + periapsis_m) / 2
        five_periods_s = 5 * 2π * sqrt(semimajor_axis_m^3 / 3.986004415e14)
        output_root = mktempdir()
        summary = run_corridor_case(config, "frequent", output_root;
            horizon_s=five_periods_s, smoke=false)

        @test summary.horizon_completed
        @test !summary.full_mission_verified
        burns = CSV.read(joinpath(_case_output(output_root, "frequent"), "burns.csv"), DataFrame)
        completed = filter(:execution_complete => ==(true), burns)
        @test nrow(completed) >= 2
        @test length(unique(completed.target_pass)) >= 2
        @test length(unique(completed.target_pass)) == nrow(completed)
        @test summary.maintenance_burn_count == nrow(completed)
        @test summary.delivered_delta_v_mps ≈ sum(completed.delivered_delta_v_mps) rtol=1e-10
        @test sum(completed.requested_delta_v_mps) <=
            Float64(config.values["comparison"]["total_delta_v_budget_mps"]) + 1e-9
        @test all(completed.propellant_used_kg .> 0.0)
        @test all(completed.mass_at_start_kg .> completed.mass_at_stop_kg)
        @test all(isapprox.(completed.propellant_used_kg,
            completed.mass_at_start_kg .- completed.mass_at_stop_kg; rtol=1e-8, atol=1e-10))

        for burn in eachrow(completed)
            rocket_equation_delta_v = 9.80665 * Float64(config.values["vehicle"]["specific_impulse_s"]) *
                log(burn.mass_at_start_kg / burn.mass_at_stop_kg)
            @test isapprox(burn.delivered_delta_v_mps, rocket_equation_delta_v; rtol=1e-6, atol=1e-6)
            @test isapprox(burn.delivered_delta_v_mps, burn.requested_delta_v_mps; rtol=1e-4, atol=0.02)
            @test burn.execution_stop_time_s > burn.execution_start_time_s
        end

        ordered = sort(collect(eachrow(completed)); by=row -> row.execution_start_time_s)
        for idx in 2:length(ordered)
            @test ordered[idx - 1].execution_stop_time_s <= ordered[idx].execution_start_time_s + 1e-6
        end
        passes = CSV.read(joinpath(_case_output(output_root, "frequent"), "passes.csv"), DataFrame)
        n_pass = nrow(passes)
        verified = filter(:target_pass => <=(n_pass), completed)
        trailing = filter(:target_pass => >(n_pass), completed)
        @test length(unique(verified.target_pass)) >= 2
        for burn in eachrow(verified)
            target_row = findfirst(==(burn.target_pass), passes.pass)
            @test target_row !== nothing
            if target_row !== nothing
                peri_m = passes.physical_periapsis_altitude_m[target_row]
                @test 95_000.0 <= peri_m <= 100_000.0
                @test isapprox(burn.achieved_periapsis_altitude_m, peri_m; atol=1e-6)
            end
        end
        # A trailing burn is allowed: the horizon is fixed in initial periods while drag shortens the period,
        # so a burn at the last apoapsis can complete before its target periapsis occurs.
        @test nrow(trailing) <= 1
        for burn in eachrow(trailing)
            @test burn.target_pass == n_pass + 1
            @test ismissing(burn.achieved_periapsis_altitude_m) || isnan(burn.achieved_periapsis_altitude_m)
            @test burn.execution_start_time_s > passes.time_s[end]
        end
    end

    @testset "safe near-passive guidance coasts without thrust or propellant loss" begin
        config = _fresh_corridor_config()
        apoapsis_m = 6_378_136.6 + Float64(config.values["orbit"]["apoapsis_altitude_m"])
        periapsis_m = 6_378_136.6 + Float64(config.values["orbit"]["initial_periapsis_altitude_m"])
        semimajor_axis_m = (apoapsis_m + periapsis_m) / 2
        two_periods_s = 2 * 2π * sqrt(semimajor_axis_m^3 / 3.986004415e14)
        output_root = mktempdir()
        summary = run_corridor_case(config, "near_passive", output_root;
            horizon_s=two_periods_s, smoke=false)

        @test summary.horizon_completed
        @test summary.maintenance_burn_count == 0
        @test summary.delivered_delta_v_mps == 0.0
        @test summary.propellant_remaining_kg == Float64(config.values["vehicle"]["propellant_mass_kg"])
        burns = CSV.read(joinpath(_case_output(output_root, "near_passive"), "burns.csv"), DataFrame)
        @test nrow(burns) == 0
        passes = CSV.read(joinpath(_case_output(output_root, "near_passive"), "passes.csv"), DataFrame)
        @test nrow(passes) >= 2
        @test any(passes.decision .== "inside_safety_band")
        @test maximum(passes.mass_kg) - minimum(passes.mass_kg) <= 1e-9
    end
end

end # module PeriapsisCorridorRuntimeTests
