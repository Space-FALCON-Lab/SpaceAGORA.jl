module HYPRExplorationScoreTests
using Test, StaticArrays, LinearAlgebra, Random
using Base.Threads: @threads, maxthreadid, threadid

const GNC = joinpath(@__DIR__, "..", "..", "src", "gnc")
include(joinpath(GNC, "hypr", "hypr_utils.jl"))
using .HYPRUtils
for file in ("reference_geometry/station_geometry.jl", "reference_geometry/cubesat_geometry.jl",
             "reference_geometry/rpo_reference_geometry.jl", "distances/mesh_distance.jl", "distances/clearance.jl")
    include(joinpath(GNC, "navigation", "rpo_nav", file))
end
for file in ("pso_parameters.jl", "path_retiming.jl", "path_sampling.jl", "path_costs.jl",
             "pso_adaptive_policy.jl", "pso_helpers.jl", "pso_refinement.jl", "rrt_connect.jl", "pso_path_planning.jl")
    include(joinpath(GNC, "guidance", "rpo", "hypr", file))
end

@testset "HYPR RRT exploration score" begin
    geometry = RPOReferenceGeometry(RPOStationGeometry(reshape([100.0, 0.0, 0.0], 3, 1)))
    start, goal = (0.0, 0.0, 0.0), (2.0, 0.0, 0.0)
    cfg = RPOPSOConfig(
        n_waypoints=1, n_particles=4, n_iters=1, curve_type=:polyline,
        adaptive_n_waypoints_min=1, adaptive_n_particles_min=4, adaptive_n_iters_min=1,
        adaptive_effort_min_fraction=1.0, adaptive_effort_max_fraction=1.0,
        rrt_warmstart_iters=100, refinement_enable=false, cull_enable=false,
        reexplore_enable=false, stagnation_learning_enable=false,
    )
    direct = hcat(collect(start), collect(goal))
    # Length 3 m versus a 2 m direct route: 50% excess length, 25 iterations used.
    detour_path = [0.0 1.0 2.0; 0.0 sqrt(1.25) 0.0; 0.0 0.0 0.0]
    diagnostics = (attempted=true, path_found=true, iterations=25)
    _, score = rpo_adaptive_pso_config(cfg, start, goal, geometry;
        rrt_path=detour_path, rrt_diagnostics=diagnostics)
    @test score.detour_ratio ≈ 1.5
    @test score.detour ≈ 0.5
    @test score.search_effort ≈ 1 - exp(-0.25)
    @test score.explore ≈ (1 + 1 - exp(-0.25)) / 3
    _, larger_budget = rpo_adaptive_pso_config(
        rpo_pso_config(cfg; rrt_warmstart_iters=25_000), start, goal, geometry;
        rrt_path=detour_path, rrt_diagnostics=diagnostics)
    @test larger_budget.search_effort == score.search_effort
    @test larger_budget.explore == score.explore
    @test RPOPSOConfig().rrt_warmstart_iters == 25_000
    @test RPOPSORRTConnectWarmstartSettings().n_iters == 25_000
    @test RPORRTConnectSettings().n_iters == 25_000
    @test RPOPSOConfig().rrt_warmstart_runtime_limit_s == Inf

    # A short route can still require many search iterations.
    _, hard_direct = rpo_adaptive_pso_config(cfg, start, goal, geometry;
        rrt_path=direct, rrt_diagnostics=merge(diagnostics, (iterations=100,)))
    @test hard_direct.explore ≈ (1 - exp(-1.0)) / 3
    @test hard_direct.search_effort > score.search_effort
    _, saturated = rpo_adaptive_pso_config(cfg, start, goal, geometry;
        rrt_path=[0.0 1.0 2.0; 0.0 10.0 0.0; 0.0 0.0 0.0],
        rrt_diagnostics=merge(diagnostics, (iterations=200,)))
    @test saturated.explore ≈ (2 + 1 - exp(-2.0)) / 3
    @test hard_direct.search_effort < saturated.search_effort < 1.0
    _, failed = rpo_adaptive_pso_config(cfg, start, goal, geometry;
        rrt_diagnostics=merge(diagnostics, (path_found=false, iterations=1,)))
    @test failed.explore == 1.0

    # With fixed geometry, Bezier control counts follow exploration and its bounds.
    bezier_cfg = rpo_pso_config(cfg; curve_type=:bezier)
    easy_cfg, _ = rpo_adaptive_pso_config(bezier_cfg, start, goal, geometry;
        rrt_path=direct, rrt_diagnostics=merge(diagnostics, (iterations=0,)))
    detour_cfg, _ = rpo_adaptive_pso_config(bezier_cfg, start, goal, geometry;
        rrt_path=detour_path, rrt_diagnostics=diagnostics)
    failed_cfg, _ = rpo_adaptive_pso_config(bezier_cfg, start, goal, geometry;
        rrt_diagnostics=merge(diagnostics, (path_found=false,)))
    @test easy_cfg.n_waypoints == 1
    @test detour_cfg.n_waypoints == 2
    @test failed_cfg.n_waypoints == 4
    capped_cfg, _ = rpo_adaptive_pso_config(
        rpo_pso_config(bezier_cfg; adaptive_n_waypoints_max=2), start, goal, geometry;
        rrt_diagnostics=merge(diagnostics, (path_found=false,)))
    @test capped_cfg.n_waypoints == 2

    # Exercise actual RRT and PSO integration, with and without warm-start seeding.
    for warmstart_enabled in (false, true)
        plan = rpo_pso_plan_path(start, goal, geometry,
            rpo_pso_config(cfg; rrt_warmstart_enable=warmstart_enabled); rng=MersenneTwister(7))
        @test plan.adaptive.explore == 0.0
        @test plan.adaptive.detour_ratio == 1.0
        @test plan.adaptive.search_effort == 0.0
        @test plan.warmstart.attempted == warmstart_enabled
        @test isfinite(plan.cost)
    end
    _, coincident = rpo_adaptive_pso_config(cfg, start, start, geometry; rng=MersenneTwister(7))
    @test coincident.explore == 0.0
end
end
