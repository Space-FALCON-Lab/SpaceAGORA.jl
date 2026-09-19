using Test
using LinearAlgebra
using Random
using StaticArrays
using CSV
using DataFrames
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel

# The example is a script with its own includes and constants; load it into a
# dedicated module so nothing leaks into Main or other test files.
module RPOStationExample
using SpaceAGORA
end
Base.include(RPOStationExample, normpath(joinpath(@__DIR__, "..", "..", "..", "examples", "Earth_RPO_CubeSat_MPC.jl")))
const RPOX = RPOStationExample
module RPOISSDemo
using SpaceAGORA
end
const RPO_STATION_SPICE = joinpath(RPOX.REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")

# A small swarm and a sparse cloud keep these tests short: they check the
# station plumbing and bounded behaviour, not planner quality.
const RPO_STATION_SMALL = (pso_n_particles=8, pso_n_iters=2, n_station_points=256, verbose=false)

# Deviation at which the replanning example (`vehicle_displacement_*`
# scenarios in Earth_RPO_CubeSat_MPC_Replanning.jl) retimes the plan. A flown
# hop that stays below it never needs the shipped retime policy.
const RPO_RETIME_DEVIATION_M = 0.75

"""Points on the six faces of a box with the given half extents, 3 x N."""
function box_shell_pointcloud(half_extents; n::Integer=7)
    hx, hy, hz = Float64.(half_extents)
    us = range(-1.0, 1.0; length=n)
    cols = SVector{3, Float64}[]
    for a in us, b in us
        push!(cols, SVector(hx, a * hy, b * hz))
        push!(cols, SVector(-hx, a * hy, b * hz))
        push!(cols, SVector(a * hx, hy, b * hz))
        push!(cols, SVector(a * hx, -hy, b * hz))
        push!(cols, SVector(a * hx, b * hy, hz))
        push!(cols, SVector(a * hx, b * hy, -hz))
    end
    return Matrix{Float64}(reduce(hcat, cols))
end

station_link(demo) = demo.args.dynamics_model.spacecraft[2].root

function build_station_demo(; kwargs...)
    settings = merge(RPO_STATION_SMALL, (results_directory=mktempdir(),), NamedTuple(kwargs))
    return RPOX.build_rpo_cubesat_mpc_demo(; settings...)
end

if !isdir(RPO_STATION_SPICE)
    @info "SPICE kernels absent; skipping the configurable-station RPO tests" RPO_STATION_SPICE
else
    @testset "RPO configurable station" begin
        @testset "defaults reproduce the Gateway scenario" begin
            d0 = build_station_demo()
            d1 = build_station_demo(;
                station_points=nothing,
                station_keepout_radius_m=0.25,
                station_name="gateway_core",
                station_dims_m=(4.0, 2.0, 2.0),
                station_mass_kg=500.0,
                station_ref_area_m2=nothing,
                safe_distance_m=0.1,
                cost_ref_distance_m=20.0,
                search_margin_m=nothing,
                sample_ds_m=0.05,
            )
            for d in (d0, d1)
                link = station_link(d)
                @test link.m == 500.0
                @test Tuple(link.dims) == (4.0, 2.0, 2.0)
                @test link.ref_area == 8.0
                @test d.geometry.station.name == "gateway_core"
                @test d.geometry.station.keepout_radius_m == 0.25
                @test size(d.geometry.station.points_body) == (3, 256)
                @test d.station.source === :gateway
                @test d.station.name == "gateway_core"
                @test d.station.n_points == 256
                @test d.station.dims_m == (4.0, 2.0, 2.0)
                @test d.station.mass_kg == 500.0
                @test d.station.ref_area_m2 == 8.0
                @test d.station.safe_distance_m == 0.1
                @test d.guidance.safe_distance_m == 0.1
                @test d.pso_config.sample_ds_m == 0.05
                @test d.pso_config.cost_ref_distance_m == 20.0
                @test d.pso_config.search_margin_m == SM.RPOPSOConfig().search_margin_m
                @test d.station.search_margin_m == d.pso_config.search_margin_m
            end
            # Same seeds through the new keywords give the same sampled cloud.
            @test d0.geometry.station.points_body == d1.geometry.station.points_body
            # Identical seeds reproduce the plan at every supported thread count.
            @test d0.plan_result.path == d1.plan_result.path
            @test d0.plan_result.cost == d1.plan_result.cost
            @test d0.initial_plan.t_ref_s == d1.initial_plan.t_ref_s
            @test d0.initial_plan.r_ref_rtn == d1.initial_plan.r_ref_rtn
            @test d0.initial_plan.v_ref_rtn == d1.initial_plan.v_ref_rtn
            capped = build_station_demo(; reference_max_speed_mps=0.05)
            @test capped.pso_config.retime_max_speed_mps == 0.05
            @test capped.plan_result.path == d0.plan_result.path
            @test maximum(norm.(eachcol(capped.initial_plan.v_ref_rtn))) <= 0.05 + 1e-10
        end

        @testset "custom point cloud, dimensions and mass" begin
            cloud = box_shell_pointcloud((1.0, 0.5, 0.5))   # a 2 x 1 x 1 m box about the origin
            start = SVector{3, Float64}(-4.0, -2.0, 1.0)
            goal = SVector{3, Float64}(-3.0, 2.0, 0.5)      # same side of the box as the start
            dc = build_station_demo(;
                station_points=cloud,
                station_keepout_radius_m=0.3,
                station_name="box",
                station_dims_m=(2.0, 1.0, 1.0),
                station_mass_kg=1200.0,
                safe_distance_m=0.2,
                cost_ref_distance_m=12.0,
                search_margin_m=6.0,
                sample_ds_m=0.1,
                start_rtn=start,
                goal_rtn=goal,
            )
            @test dc.station.source === :custom
            @test dc.station.name == "box"
            @test dc.station.n_points == size(cloud, 2)
            @test dc.geometry.station.points_body == cloud
            @test dc.geometry.station.name == "box"
            @test dc.geometry.station.keepout_radius_m == 0.3
            link = station_link(dc)
            @test link.m == 1200.0
            @test Tuple(link.dims) == (2.0, 1.0, 1.0)
            @test link.ref_area == 1.0                       # y-z face of the custom box
            @test dc.station.ref_area_m2 == 1.0
            @test link.inertia ≈ (1200.0 / 12) * Diagonal([1.0 + 1.0, 4.0 + 1.0, 4.0 + 1.0])
            @test dc.station.dims_m == (2.0, 1.0, 1.0)
            @test dc.station.mass_kg == 1200.0
            @test dc.pso_config.sample_ds_m == 0.1
            @test dc.pso_config.cost_ref_distance_m == 12.0
            @test dc.pso_config.search_margin_m == 6.0
            @test dc.guidance.safe_distance_m == 0.2
            @test dc.station.safe_distance_m == 0.2
            @test dc.station.cost_ref_distance_m == 12.0
            @test dc.station.search_margin_m == 6.0
            @test dc.station.sample_ds_m == 0.1

            # An explicit area wins over the derived face area.
            de = build_station_demo(; station_points=cloud, station_dims_m=(2.0, 1.0, 1.0), station_ref_area_m2=7.5)
            @test station_link(de).ref_area == 7.5
            @test de.station.ref_area_m2 == 7.5

            # Bounded planner: the plan runs from start to goal and both the
            # sampled path and the retimed reference clear the keep-out surface
            # of the custom cloud by the requested safe distance.
            @test dc.initial_plan.r_ref_rtn[:, 1] ≈ start atol=1e-9
            @test dc.initial_plan.r_ref_rtn[:, end] ≈ goal atol=1e-9
            sampled = SM.GuidanceHooks.rpo_sample_path(dc.initial_plan.path_rtn, dc.pso_config.sample_ds_m; curve_type=dc.pso_config.curve_type)
            for samples in (sampled, dc.initial_plan.r_ref_rtn)
                stats = SM.rpo_path_clearance_stats(samples, dc.geometry; safe_distance_m=0.2)
                @test stats.violation_count == 0
                @test stats.min_clearance >= 0.2
            end
        end

        @testset "invalid station inputs are refused" begin
            cloud = box_shell_pointcloud((1.0, 0.5, 0.5))
            @test_throws ArgumentError build_station_demo(; station_points=cloud[1:2, :])
            @test_throws ArgumentError build_station_demo(; station_dims_m=(2.0, 0.0, 1.0))
            @test_throws ArgumentError build_station_demo(; station_dims_m=(2.0, 1.0))
            @test_throws ArgumentError build_station_demo(; station_mass_kg=-1.0)
            @test_throws ArgumentError build_station_demo(; station_ref_area_m2=0.0)
            @test_throws ArgumentError build_station_demo(; station_keepout_radius_m=-0.1)
            @test_throws ArgumentError build_station_demo(; reference_max_speed_mps=0.0)
            @test_throws ArgumentError build_station_demo(; reference_max_speed_mps=Inf)
        end

        @testset "bounded tracking of a short hop beside the default station" begin
            start = SVector{3, Float64}(-8.0, -4.0, 2.0)
            goal = SVector{3, Float64}(-5.5, -2.5, 1.0)
            dt = build_station_demo(; mission_time=30.0, start_rtn=start, goal_rtn=goal, data_rate_s=1.0)
            run_simulation(dt.args)
            csv = joinpath(dt.args.simulation_settings.results_directory, "simulation_results.csv")
            @test isfile(csv)
            df, actual_rtn, ref_rtn, err = RPOX._rpo_postprocess(csv, dt)
            @test nrow(df) >= 10
            @test all(isfinite, err)
            @test df.time[end] >= dt.initial_plan.t_ref_s[end]
            # Bounded tracking: the flown deviation from the reference stays
            # below the retime threshold of the shipped replanning example and
            # below the clearance the plan keeps from the keep-out surface, so
            # the chaser never crosses that surface, and it ends on the goal
            # within the planning safe distance.
            planned = SM.rpo_path_clearance_stats(dt.initial_plan.r_ref_rtn, dt.geometry; safe_distance_m=0.0)
            @test maximum(err) < RPO_RETIME_DEVIATION_M
            @test maximum(err) < planned.min_clearance
            @test SM.rpo_path_clearance_stats(actual_rtn, dt.geometry; safe_distance_m=0.0).violation_count == 0
            @test norm(actual_rtn[:, end] - goal) < dt.station.safe_distance_m
            @test df.sc1_mass[end] < df.sc1_mass[1]     # propellant was spent
        end

        @testset "ISS demo provenance and reuse rules" begin
            # Load the demo script without running it (main() only runs as a program).
            demo_out = mktempdir()
            withenv("SPACEAGORA_VIEWER_DEMO_OUT" => demo_out, "SPACEAGORA_DEMO_SMOKE" => "0") do
                Base.include(RPOISSDemo, normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "dev", "viewer_demos", "iss_hypr.jl")))
            end
            D = RPOISSDemo
            full = D.iss_hypr_inputs(; smoke=false)
            smoke = D.iss_hypr_inputs(; smoke=true)
            @test full.smoke === false && smoke.smoke === true
            @test full.station_dims_m == (73.0, 109.0, 20.0)
            @test full.reference_max_speed_mps == smoke.reference_max_speed_mps == 0.1
            @test length(full.model_sha256) == 64
            @test D.iss_hypr_outdir(full) == D.iss_hypr_outdir(D.iss_hypr_inputs(; smoke=false))   # same inputs, same directory
            @test D.iss_hypr_outdir(full) != D.iss_hypr_outdir(smoke)                                # different inputs, different directory
            @test startswith(D.iss_hypr_outdir(full), demo_out)

            outdir = D.iss_hypr_outdir(smoke)
            mkpath(outdir)
            prefix = joinpath(outdir, "simulation_results")
            sidecar = joinpath(outdir, "iss_hypr_provenance.json")
            planfile = joinpath(outdir, "iss_hypr_plan.json")
            @test !D.iss_hypr_matching_run(sidecar, smoke, prefix, planfile)
            write(sidecar, D.JSON.json(Dict("inputs" => D._json_roundtrip(smoke))))
            @test !D.iss_hypr_matching_run(sidecar, smoke, prefix, planfile)       # scene and plan still missing
            write(prefix * "_scene.json", "{}")
            @test !D.iss_hypr_matching_run(sidecar, smoke, prefix, planfile)       # plan still missing
            write(planfile, D.JSON.json(Dict("path_rtn" => [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])))
            @test D.iss_hypr_matching_run(sidecar, smoke, prefix, planfile)
            @test !D.iss_hypr_matching_run(sidecar, full, prefix, planfile)        # other inputs never reuse this run
            @test D._matrix3(D.JSON.parsefile(planfile)["path_rtn"]) == [1.0 4.0; 2.0 5.0; 3.0 6.0]
        end
    end
end
