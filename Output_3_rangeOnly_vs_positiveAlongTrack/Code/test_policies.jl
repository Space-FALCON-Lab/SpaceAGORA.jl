import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "9_SpaceAGORA.jl-main_full integration")))

using CSV, DataFrames, Dates, LinearAlgebra, Random, StaticArrays, Test, TOML
include("simulation.jl")
using .ComparisonSimulation
using .ComparisonSimulation.ComparisonPolicies

# Host the unchanged reference module without loading its unrelated simulation engine.
module ReferenceSpaceAGORA
    module AbstractTypes
        abstract type AbstractForceTorqueModel end
    end
    module Dynamics
        module DynamicEffectors
            function calcForceTorque end
            function solver_partition end
        end
        const SOURCE = normpath(joinpath(@__DIR__, "..", "..", "..", "ORACLE_SpaceAGORA_Integration",
            "2_SpaceAGORA.jl", "src", "dynamics", "coupled", "force_torque_models",
            "laser_link_effectors.jl"))
        include(SOURCE)
    end
end
const Reference = ReferenceSpaceAGORA.Dynamics.LaserLinkEffectors
const CS = ComparisonSimulation

function change_settings(s::Settings; kwargs...)
    existing = (; (name => getfield(s, name) for name in fieldnames(Settings))...)
    Settings(; merge(existing, (; kwargs...))...)
end

const VALIDATION = VALIDATION_ROOT
mkpath(VALIDATION)
convergence = NamedTuple[]

@testset "Range-only and positive along-track comparison" begin
    @testset "Isolated code and paper output layout" begin
        @test basename(CODE_ROOT) == "Code"
        @test DEFAULT_OUTPUT == normpath(joinpath(CODE_ROOT, "..", "Output"))
        @test dirname(CS.GRAVITY_SOURCE) == CODE_ROOT
        @test VALIDATION == normpath(joinpath(CODE_ROOT, "..", "..", "Tests",
            "Output_3_rangeOnly_vs_positiveAlongTrack", "validation"))
        for policy in (:naive_next_entering, :positive_along_track)
            @test relpath(case_directory(DEFAULT_OUTPUT, Settings(), policy), DEFAULT_OUTPUT) ==
                joinpath("paper_plot_mode", "h1050km_t1000km", "ih0.0deg_it0.5deg",
                    "N200", "T504570s", "$(policy)_e0.0000_nu1.0000")
        end
    end

    @testset "Policy transitions and reference equivalence" begin
        target = SVector(7_000_000.0, 0.0, 0.0)
        velocities = fill(SVector(0.0, 7500.0, 0.0), 4)
        positions(offsets) = [target, (target + SVector{3,Float64}(o) for o in offsets)...]
        range_policy = LaserPolicy(:naive_next_entering; helper_indices=[2, 3, 4], range_m=100)
        @test update_schedule!(range_policy,
            positions([(40, 0, 0), (50, 0, 0), (200, 0, 0)]), velocities) == 2
        @test update_schedule!(range_policy,
            positions([(80, 0, 0), (10, 0, 0), (200, 0, 0)]), velocities) == 2
        @test update_schedule!(range_policy,
            positions([(150, 0, 0), (10, 0, 0), (200, 0, 0)]), velocities) == 0
        @test update_schedule!(range_policy,
            positions([(150, 0, 0), (10, 0, 0), (200, 0, 0)]), velocities) == 0
        @test update_schedule!(range_policy,
            positions([(150, 0, 0), (10, 0, 0), (100, 0, 0)]), velocities) == 4
        @test range_policy.activations == 2

        along = LaserPolicy(:positive_along_track; helper_indices=[2, 3, 4], range_m=100)
        @test update_schedule!(along,
            positions([(80, -60, 0), (0, -50, 0), (200, 0, 0)]), velocities) == 3
        @test update_schedule!(along,
            positions([(0, -50, 0), (80, -10, 0), (200, 0, 0)]), velocities) == 3
        @test update_schedule!(along,
            positions([(0, -50, 0), (0, 50, 0), (200, 0, 0)]), velocities) == 2
        @test update_schedule!(along,
            positions([(0, 50, 0), (0, 50, 0), (80, 0, 0)]), velocities) == 0
        @test update_schedule!(along,
            positions([(0, -100, 0), (0, -100, 0), (0, 0, 0)]), velocities) == 2

        for schedule in (:naive_next_entering, :positive_along_track)
            rebuilt = LaserPolicy(schedule; helper_indices=[2, 3, 4], range_m=100)
            original = Reference.OpenCavityLaserLinkModel(; schedule,
                helper_indices=[2, 3, 4], range_m=100)
            rng = MersenneTwister(2027)
            for step in 1:1000
                offsets = [randn(rng, SVector{3,Float64}) * 65 for _ in 1:3]
                step % 17 == 0 && (offsets[1] = SVector(0.0, -100.0, 0.0))
                step % 19 == 0 && (offsets[2] = offsets[1])
                pos = positions(offsets)
                @test update_schedule!(rebuilt, pos, velocities) ==
                    Reference.update_laser_link_schedule!(original, pos, velocities)
                @test rebuilt.previous_in_range == original.previous_in_range
                @test rebuilt.activations == original.link_activation_count
                @test pair_force(rebuilt, pos[1], pos[2]) ≈
                    Reference.laser_link_pair_force(original, pos[1], pos[2])
            end
        end
        @test iszero(pair_force(along, target, target))
        @test iszero(pair_force(along, target, target + SVector(101.0, 0.0, 0.0)))
        @test norm(pair_force(along, target, target + SVector(100.0, 0.0, 0.0))) ≈ along.force_N
        @test_throws ArgumentError LaserPolicy(:gve_sma; helper_indices=[2])
        @test_throws ArgumentError LaserPolicy(:naive_next_entering; helper_indices=[1])
        @test_throws ArgumentError LaserPolicy(:naive_next_entering; helper_indices=[2, 2])
        @test_throws ArgumentError LaserPolicy(:positive_along_track; helper_indices=[2], range_m=NaN)
    end

    @testset "Requested geometry and physical forces" begin
        s = Settings()
        state = initial_state(s)
        @test size(state) == (6, 201)
        @test s.orbits == 80
        @test s.j2 == J2 > 0
        @test s.target_inclination_deg - s.helper_inclination_deg == 0.5
        @test duration(s) == 80target_period(s)
        for index in axes(state, 2)
            target = index == 1
            r, v = CS.position(state, index), CS.velocity(state, index)
            expected_radius = EARTH_RADIUS_M + (target ? 1_000_000 : 1_050_000)
            @test norm(r) ≈ expected_radius atol=1e-8 rtol=0
            @test norm(v) ≈ sqrt(MU / expected_radius) atol=1e-10 rtol=0
            oe = orbital_elements(r, v)
            @test oe.eccentricity < 1e-14
            @test oe.inclination_deg ≈ (target ? 0.5 : 0.0) atol=1e-12
        end
        @test atan(state[2, 2], state[1, 2]) == 0
        @test atan(state[2, 3], state[1, 3]) ≈ deg2rad(360 / 200)
        @test orbital_elements(CS.position(state, 1), CS.velocity(state, 1)).raan_deg ≈ 0 atol=1e-12

        planet = (μ=MU, J2=J2, Rp_e=EARTH_RADIUS_M)
        radius = EARTH_RADIUS_M + 1_000_000
        equatorial = CS._inverse_squared_j2_gravity_accel(SVector(radius, 0.0, 0.0), planet)
        polar = CS._inverse_squared_j2_gravity_accel(SVector(0.0, 0.0, radius), planet)
        @test equatorial[1] ≈ -MU / radius^2 * (1 + 1.5J2 * (EARTH_RADIUS_M / radius)^2)
        @test polar[3] ≈ -MU / radius^2 * (1 - 3J2 * (EARTH_RADIUS_M / radius)^2)

        policy = LaserPolicy(:naive_next_entering; helper_indices=collect(2:201))
        params = CS.PropagatorParameters(s, policy, planet)
        baseline, with_laser = similar(state), similar(state)
        CS.dynamics!(baseline, state, params, 0.0)
        policy.active_helper_idx = 2
        CS.dynamics!(with_laser, state, params, 0.0)
        laser = (with_laser[4:6, :] - baseline[4:6, :]) * s.mass_kg
        @test norm(laser[:, 1]) ≈ CS.force_magnitude(s) atol=1e-12
        @test norm(sum(laser; dims=2)) < 1e-12
        @test laser[:, 1] ≈ -laser[:, 2] atol=1e-12
        @test all(iszero, laser[:, 3:end])
    end

    @testset "Conservation, J2 precession and step refinement" begin
        s = Settings(helpers=2, orbits=2.0, power_W=0.0, j2=0.0, samples_per_orbit=100)
        central = propagate(s, :naive_next_entering)
        elements = target_table(central)
        @test maximum(abs, elements.delta_a_km) < 1e-5
        @test maximum(elements.eccentricity) < 1e-9
        @test maximum(abs, elements.delta_inclination_deg) < 1e-10
        @test maximum(abs, elements.delta_raan_deg) < 1e-9

        j2_result = propagate(change_settings(s; j2=J2), :naive_next_entering)
        j2_table = target_table(j2_result)
        n = 2pi / target_period(s)
        a = EARTH_RADIUS_M + s.target_altitude_km * 1000
        predicted_raan = rad2deg(-1.5J2 * n * (EARTH_RADIUS_M / a)^2 *
            cosd(s.target_inclination_deg) * duration(s))
        @test j2_table.delta_raan_deg[end] < 0
        @test abs(j2_table.delta_raan_deg[end] - predicted_raan) < 0.01
        energies = map(j2_result.solution.u) do u
            r, v = CS.position(u, 1), CS.velocity(u, 1)
            radius = norm(r)
            potential = -MU / radius * (1 - J2 * (EARTH_RADIUS_M / radius)^2 *
                (3(r[3] / radius)^2 - 1) / 2)
            dot(v, v) / 2 + potential
        end
        @test maximum(abs, energies .- first(energies)) / abs(first(energies)) < 1e-10

        s = Settings(orbits=2.0, samples_per_orbit=100, max_step_s=10.0)
        for policy in (:naive_next_entering, :positive_along_track)
            coarse = propagate(s, policy)
            fine = propagate(change_settings(s; max_step_s=5.0), policy)
            coarse_table, fine_table = target_table(coarse), target_table(fine)
            differences = (max_a_difference_m=1000maximum(abs, coarse_table.a_km - fine_table.a_km),
                max_e_difference=maximum(abs, coarse_table.eccentricity - fine_table.eccentricity),
                max_i_difference_deg=maximum(abs, coarse_table.inclination_deg - fine_table.inclination_deg),
                max_raan_difference_deg=maximum(abs, coarse_table.raan_deg - fine_table.raan_deg))
            push!(convergence, merge((policy=string(policy), orbits=2.0,
                coarse_max_step_s=10.0, fine_max_step_s=5.0), differences))
            @test coarse.policy.activations > 0
            @test fine.policy.activations > 0
            @test differences.max_a_difference_m < 5.0
            @test differences.max_e_difference < 1e-6
            @test differences.max_i_difference_deg < 1e-5
            @test differences.max_raan_difference_deg < 1e-3
        end
    end

    @testset "Persistent output and overwrite protection" begin
        mktempdir() do root
            s = Settings(helpers=3, orbits=0.02)
            infos = run_comparison(; settings=s, root)
            @test all(info -> info["status"] == "verified", infos)
            @test all(info -> info["trajectory_rows"] == 5, infos)
            @test all(info -> info["settings"]["j2"] == J2, infos)
            for policy in (:naive_next_entering, :positive_along_track)
                directory = case_directory(root, s, policy)
                @test isfile(joinpath(directory, "trajectory.feather"))
                @test isfile(joinpath(directory, "target_orbital_elements.csv"))
                @test isfile(joinpath(directory, "schedule.csv"))
                @test isfile(joinpath(directory, "run_info.toml"))
            end
            @test_throws ErrorException run_comparison(; settings=s, root)
        end
    end
end

CSV.write(joinpath(VALIDATION, "two_orbit_step_refinement.csv"), DataFrame(convergence))
open(joinpath(VALIDATION, "tests.toml"), "w") do io
    TOML.print(io, Dict("passed" => true, "finished_utc" => string(now(UTC)),
        "reference_policy_steps_per_policy" => 1000,
        "reference_source" => ReferenceSpaceAGORA.Dynamics.SOURCE,
        "reference_source_sha256" => CS.filehash(ReferenceSpaceAGORA.Dynamics.SOURCE),
        "production_settings" => settings_dict(Settings()),
        "step_refinement_scope" => "Two orbits, 200 helpers, J2; 10 s versus 5 s maximum steps"))
end
