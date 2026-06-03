using Test
using StaticArrays
using LinearAlgebra
using Random
using SpaceAGORA
const SM = SpaceAGORA.SimulationModel

@testset "RPO RTN frame helpers" begin
    r_t = SVector{3, Float64}(7000e3, 0.0, 0.0)
    v_t = SVector{3, Float64}(0.0, 7500.0, 0.0)
    r_rel = SVector{3, Float64}(10.0, -2.0, 1.0)
    v_rel = SVector{3, Float64}(0.01, 0.02, -0.01)
    r_c, v_c = SM.GuidanceHooks.rtn_to_inertial_relative_state(r_rel, v_rel, r_t, v_t)
    x = SM.GuidanceHooks.inertial_to_rtn_relative_state(r_c, v_c, r_t, v_t)
    @test x[1:3] ≈ r_rel atol=1e-9
    @test x[4:6] ≈ v_rel atol=1e-9
end

@testset "RPO geometry and PSO basics" begin
    @test isfile(SpaceAGORA.station_cad_path(:gateway))
    gateway_points = SpaceAGORA.load_rpo_station_cad_pointcloud(:gateway; n_points=16, rng=MersenneTwister(740))
    @test size(gateway_points) == (3, 16)
    points = SpaceAGORA.load_rpo_station_pointcloud(:demo)
    station = SM.RPOStationGeometry(points; keepout_radius_m=0.25)
    geom = SM.RPOReferenceGeometry(station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.2, 0.2, 0.3)))
    stats = SM.rpo_clearance_to_station(SVector{3, Float64}(8.0, 0.0, 0.0), geom)
    @test stats.clearance > 0.0
    @test SM.GuidanceHooks.rpo_obstacle_sigmoid_penalty(0.08, 0.09, 1.0e6) > 1.0 - 1.0e-9
    @test SM.GuidanceHooks.rpo_obstacle_sigmoid_penalty(0.10, 0.09, 1.0e6) < 1.0e-9

    sigmoid_geom = SM.RPOReferenceGeometry(
        SM.RPOStationGeometry(reshape([0.0, 0.0, 0.0], 3, 1); keepout_radius_m=0.0);
        chaser=SM.RPOCubeSatGeometry(dims_m=(0.02, 0.02, 0.02)),
    )
    sigmoid_samples = [
        0.08 0.12;
        0.0  0.0;
        0.0  0.0
    ]
    sigmoid_stats = SM.GuidanceHooks.rpo_clearance_stats_from_samples(
        sigmoid_samples,
        sigmoid_geom,
        0.10;
        obstacle_sigmoid_k=1.0e6,
        obstacle_sigmoid_tol_m=0.01,
    )
    @test sigmoid_stats.violation_count == 0
    @test sigmoid_stats.obstacle_score ≈ 1.0 atol=1.0e-9
    intersect_stats = SM.GuidanceHooks.rpo_clearance_stats_from_samples(
        reshape([0.005, 0.0, 0.0], 3, 1),
        sigmoid_geom,
        0.10;
        obstacle_sigmoid_k=1.0e6,
        obstacle_sigmoid_tol_m=0.01,
    )
    @test intersect_stats.violation_count == 1
    simple_station = SM.RPOStationGeometry([0.0 2.0 -3.0; 0.0 0.0 4.0; 0.0 0.0 0.0])
    nearest, distance, idx = SM.nearest_station_point(SVector{3, Float64}(1.6, 0.2, 0.0), simple_station)
    @test idx == 2
    @test nearest ≈ SVector{3, Float64}(2.0, 0.0, 0.0)
    @test distance ≈ sqrt(0.2)

    cfg = SM.RPOPSOConfig(
        n_waypoints=1,
        n_particles=8,
        n_iters=2,
        sample_ds_m=0.5,
        retime_a_max_mps2=0.1,
        tf_s=20.0,
        cost_ref_distance_m=10.0,
        curve_type=:bezier,
        adaptive_enable=false,
    )
    plan = SM.rpo_pso_plan_path(
        SVector{3, Float64}(-6.0, 0.0, 0.0),
        SVector{3, Float64}(6.0, 0.0, 0.0),
        geom,
        cfg;
        safe_distance_m=0.1,
        rng=MersenneTwister(1),
    )
    @test size(plan.path, 1) == 3
    @test isfinite(plan.cost)
    @test plan.config.sample_ds_m == 0.1

    synced_cfg = SM.rpo_pso_config(SM.RPOPSOConfig(sample_ds_m=0.5, safe_distance_m=0.2))
    @test synced_cfg.sample_ds_m == synced_cfg.safe_distance_m
    @test SM.GuidanceHooks.rpo_hypr_sampling_density_m(SM.RPOPSOConfig(sample_ds_m=0.5), 0.2) == 0.2
    @test SM.GuidanceHooks.rpo_740_mpc_final_pso_config(safe_distance_m=0.2).sample_ds_m == 0.2

    adaptive_sampling_cfg = SM.RPOPSOConfig(
        sample_ds_m=0.05,
        adaptive_sampling_enable=true,
        adaptive_sampling_max_ds_m=0.5,
        adaptive_sampling_far_clearance_m=1.0,
        adaptive_sampling_obstacle_guard_fraction=0.5,
        curve_type=:polyline,
    )
    crossing_geom = SM.RPOReferenceGeometry(
        SM.RPOStationGeometry(reshape([0.0, 0.0, 0.0], 3, 1); keepout_radius_m=0.2);
        chaser=SM.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.1)),
    )
    crossing_path = [
        -2.0 2.0;
         0.0 0.0;
         0.0 0.0
    ]
    uniform_crossing = SM.GuidanceHooks.rpo_sample_path(crossing_path, adaptive_sampling_cfg.sample_ds_m; curve_type=:polyline)
    adaptive_crossing = SM.GuidanceHooks.rpo_sample_path(
        crossing_path,
        adaptive_sampling_cfg,
        crossing_geom;
        safe_distance_m=0.0,
        curve_type=:polyline,
    )
    @test size(adaptive_crossing, 2) < size(uniform_crossing, 2)
    @test !SM.GuidanceHooks.rpo_rrt_segment_is_safe(
        SVector{3, Float64}(-1.0, 0.0, 0.0),
        SVector{3, Float64}(1.0, 0.0, 0.0),
        crossing_geom;
        safe_distance_m=0.0,
        sample_ds_m=2.0,
        adaptive_enable=true,
    )

    retime_geom = SM.RPOReferenceGeometry(
        SM.RPOStationGeometry(reshape([0.0, 0.0, 0.0], 3, 1); keepout_radius_m=0.0);
        chaser=SM.RPOCubeSatGeometry(dims_m=(0.02, 0.02, 0.02)),
    )
    retime_cfg = SM.RPOPSOConfig(
        sample_ds_m=0.05,
        curve_type=:polyline,
        retime_dt_s=1.0,
        retime_reaction_time_s=0.0,
        retime_a_max_mps2=0.5,
        retime_speed_scale=1.0,
        retime_min_speed_mps=0.0,
        retime_max_speed_mps=Inf,
        retime_max_steps=100,
        adaptive_enable=false,
    )
    retime_path = [
        0.21 0.31;
        0.0  0.0;
        0.0  0.0
    ]
    _, _, v_zero_safe = SM.GuidanceHooks.rpo_retime_path(retime_path, retime_geom, retime_cfg; safe_distance_m=0.0)
    _, _, v_large_safe = SM.GuidanceHooks.rpo_retime_path(retime_path, retime_geom, retime_cfg; safe_distance_m=0.19)
    @test v_large_safe[1] ≈ v_zero_safe[1]
    @test v_zero_safe[1] ≈ sqrt(2.0 * retime_cfg.retime_a_max_mps2 * 0.21)

    modular_cfg = SM.rpo_pso_config(SM.RPOPSOConfigurator(
        swarm=SM.RPOPSOSwarmSettings(
            n_waypoints=2,
            n_particles=9,
            n_iters=3,
            sample_ds_m=0.5,
            curve_type=:bezier,
        ),
        objective=SM.RPOPSOObjectiveSettings(
            w_len=1.2,
            w_obs=1.0e5,
            cost_ref_distance_m=10.0,
            tf_s=20.0,
        ),
        adaptive=SM.RPOPSOAdaptiveSettings(
            enabled=false,
            n_particles_min=6,
            n_particles_max=12,
            n_iters_min=2,
            n_iters_max=5,
        ),
        cull=SM.RPOPSOCullSettings(
            enabled=true,
            fraction_max=0.25,
            start_iter=1,
        ),
        schedule=SM.RPOPSOScheduleSettings(enabled=true),
        early_stopping=SM.RPOPSOEarlyStoppingSettings(
            enabled=true,
            patience=1,
            min_iters=1,
            min_abs_improvement=1.0e99,
            require_feasible=false,
        ),
        rrt_warmstart=SM.RPOPSORRTConnectWarmstartSettings(
            enabled=true,
            n_iters=4,
        ),
        refinement=SM.RPOPSORefinementSettings(enabled=false),
        retiming=SM.RPOPSORetimingSettings(a_max_mps2=0.1),
    ))
    @test modular_cfg.n_waypoints == 2
    @test modular_cfg.n_particles == 9
    @test modular_cfg.adaptive_enable === false
    @test modular_cfg.cull_enable === true
    @test modular_cfg.schedule_enable === true
    @test modular_cfg.early_stopping_enable === true
    @test modular_cfg.early_stopping_patience == 1
    @test modular_cfg.refinement_enable === false
    @test modular_cfg.rrt_warmstart_enable === true
    @test modular_cfg.rrt_warmstart_iters == 4
    @test modular_cfg.rrt_warmstart_box_margin_m == 0.75
    @test SM.RPOPSOConfig().search_margin_m == 10.0
    @test SM.RPOPSOConfig().reexplore_search_margin_scale == 2.0
    @test SM.RPOPSOConfig().rrt_warmstart_enable === false
    warmstart_bounds_cfg = SM.rpo_pso_config(modular_cfg; pso_rrt_warmstart_box_margin=0.25)
    warmstart_lo, warmstart_hi = SM.GuidanceHooks.rpo_pso_warmstart_bounds(
        [
            -1.0 0.0 2.0;
             0.0 1.0 0.5;
             0.0 0.0 0.0
        ],
        warmstart_bounds_cfg,
    )
    @test warmstart_lo ≈ SVector{3, Float64}(-1.25, -0.25, -0.25)
    @test warmstart_hi ≈ SVector{3, Float64}(2.25, 1.25, 0.25)
    updated_cfg = SM.rpo_pso_config(modular_cfg; pso_particles=11, pso_cull_start_iter=2, pso_rrt_warmstart_iters=5)
    @test updated_cfg.n_particles == 11
    @test updated_cfg.cull_start_iter == 2
    @test updated_cfg.rrt_warmstart_iters == 5

    warm_plan = SM.rpo_pso_plan_path(
        SVector{3, Float64}(-4.0, 3.0, 0.0),
        SVector{3, Float64}(4.0, 3.0, 0.0),
        geom,
        updated_cfg;
        safe_distance_m=0.1,
        rng=MersenneTwister(2),
    )
    @test warm_plan.warmstart.enabled === true
    @test warm_plan.warmstart.attempted === true
    @test warm_plan.warmstart.path_found === true
    @test warm_plan.warmstart.n_points >= 2
    @test warm_plan.early_stopped === true
    @test length(warm_plan.cost_history) < updated_cfg.n_iters

    fixed_cfg, adaptive_stats = SM.GuidanceHooks.rpo_adaptive_pso_config(
        modular_cfg,
        SVector{3, Float64}(-6.0, 0.0, 0.0),
        SVector{3, Float64}(6.0, 0.0, 0.0),
        geom;
        safe_distance_m=0.1,
    )
    @test fixed_cfg.n_particles == modular_cfg.n_particles
    @test adaptive_stats.enabled === false

    cull_cfg = SM.RPOPSOConfig(
        n_waypoints=3,
        n_particles=4,
        n_iters=1,
        cull_enable=true,
        cull_fraction_max=0.25,
        cull_start_iter=1,
        cull_noise_scale=0.0,
        cull_arc_velocity_scale=0.0,
        adaptive_enable=false,
    )
    start = SVector{3, Float64}(0.0, 0.0, 0.0)
    goal = SVector{3, Float64}(10.0, 0.0, 0.0)
    gbest_path = [
        0.0 2.5 5.0 7.5 10.0;
        0.0 2.0 3.0 2.0  0.0;
        0.0 0.0 0.0 0.0  0.0
    ]
    gbest = [gbest_path[:, j] for j in 2:4]
    gbest = reduce(vcat, gbest)
    dim = length(gbest)
    positions = repeat(gbest, 1, 4)
    velocities = fill(1.0, dim, 4)
    pbest = copy(positions)
    pbest_cost = [1.0, 2.0, 3.0, 4.0]
    lo_rep = fill(-20.0, dim)
    hi_rep = fill(20.0, dim)
    replaced = SM.GuidanceHooks.rpo_pso_cull_swarm!(
        positions,
        velocities,
        pbest,
        pbest_cost,
        lo_rep,
        hi_rep,
        gbest,
        start,
        goal,
        3,
        cull_cfg,
        1,
        MersenneTwister(3),
    )
    @test replaced == 1
    expected_cull = zeros(dim)
    for j in 1:3
        q_star = SVector{3, Float64}(gbest_path[:, j + 1])
        q_line = (1.0 - j / 4) * start + (j / 4) * goal
        q_loc = SM.GuidanceHooks.rpo_pso_project_to_segment(q_star, gbest_path[:, j], gbest_path[:, j + 2])
        expected_cull[(3 * (j - 1) + 1):(3 * j)] .= q_star + 0.6 * (q_line - q_star) + 0.4 * (q_loc - q_star)
    end
    @test positions[:, 4] ≈ expected_cull
    @test all(velocities[:, 4] .== 0.0)

    refine_cfg = SM.RPOPSOConfig(
        n_waypoints=2,
        n_particles=4,
        n_iters=1,
        sample_ds_m=0.1,
        curve_type=:bezier,
        adaptive_enable=false,
        refinement_enable=true,
        refinement_rounds=4,
        refinement_waypoint_passes=4,
        refinement_sample_ds_m=0.05,
        retime_a_max_mps2=0.1,
        tf_s=20.0,
        cost_ref_distance_m=2.0,
    )
    far_geom = SM.RPOReferenceGeometry(
        SM.RPOStationGeometry(reshape([0.0, 10.0, 0.0], 3, 1); keepout_radius_m=0.0);
        chaser=SM.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.1)),
    )
    curved = [
        -1.0 -0.5 0.5 1.0;
         0.0  1.0 1.0 0.0;
         0.0  0.0 0.0 0.0
    ]
    before = SM.GuidanceHooks.rpo_normalized_path_cost_components(curved, far_geom, refine_cfg; safe_distance_m=0.0)
    refined, refined_cost, improved = SM.GuidanceHooks.rpo_post_refine_path(curved, far_geom, refine_cfg; safe_distance_m=0.0)
    after = SM.GuidanceHooks.rpo_normalized_path_cost_components(refined, far_geom, refine_cfg; safe_distance_m=0.0)
    @test improved === true
    @test refined_cost < before.total
    @test after.J_obs <= before.J_obs
    @test size(refined, 2) <= size(curved, 2)
end

@testset "RPO HyPR replanning decisions" begin
    station = SM.RPOStationGeometry(reshape([0.0, 10.0, 0.0], 3, 1); keepout_radius_m=0.0)
    geom = SM.RPOReferenceGeometry(station; chaser=SM.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.1)))
    path = [
        -2.0 -1.0 0.0 1.0 2.0;
         0.0  0.0 0.0 0.0 0.0;
         0.0  0.0 0.0 0.0 0.0
    ]
    plan = SM.RPOPlan(
        valid=true,
        t_ref_s=collect(0.0:1.0:4.0),
        r_ref_rtn=path,
        v_ref_rtn=zeros(3, 5),
        path_rtn=path,
        cost=0.0,
    )

    blocking = SM.GuidanceHooks.RPOReplanningConfig(
        spheres=[SM.GuidanceHooks.RPOReplanningSphere(SVector{3, Float64}(0.0, 0.0, 0.0), 0.2)],
        safe_distance_m=0.1,
        retime_clearance_m=0.5,
        sphere_surface_samples=160,
    )
    decision = SM.GuidanceHooks.rpo_replanning_decision(plan, SVector{3, Float64}(-1.5, 0.0, 0.0), geom, blocking, 1.0)
    @test decision.action == :replan

    near = SM.GuidanceHooks.RPOReplanningConfig(
        spheres=[SM.GuidanceHooks.RPOReplanningSphere(SVector{3, Float64}(0.0, 0.35, 0.0), 0.1)],
        safe_distance_m=0.1,
        retime_clearance_m=0.5,
        sphere_surface_samples=160,
    )
    decision = SM.GuidanceHooks.rpo_replanning_decision(plan, SVector{3, Float64}(-1.5, 0.0, 0.0), geom, near, 1.0)
    @test decision.action == :retime

    behind = SM.GuidanceHooks.RPOReplanningConfig(
        spheres=[SM.GuidanceHooks.RPOReplanningSphere(SVector{3, Float64}(-1.5, 0.0, 0.0), 0.3)],
        safe_distance_m=0.1,
        retime_clearance_m=0.5,
        sphere_surface_samples=160,
    )
    decision = SM.GuidanceHooks.rpo_replanning_decision(plan, SVector{3, Float64}(1.0, 0.0, 0.0), geom, behind, 1.0)
    @test decision.action == :none

    disturbed = SM.GuidanceHooks.RPOReplanningConfig(
        safe_distance_m=0.1,
        retime_clearance_m=0.5,
        tracking_error_retime_m=0.5,
    )
    decision = SM.GuidanceHooks.rpo_replanning_decision(plan, SVector{3, Float64}(0.0, 1.0, 0.0), geom, disturbed, 1.0)
    @test decision.action == :retime
    @test decision.reason == :tracking_error

    goal_changed = SM.GuidanceHooks.RPOReplanningConfig(
        desired_goal_rtn=SVector{3, Float64}(3.0, 0.5, 0.0),
        goal_change_tolerance_m=0.05,
        safe_distance_m=0.1,
        retime_clearance_m=0.5,
    )
    decision = SM.GuidanceHooks.rpo_replanning_decision(plan, SVector{3, Float64}(0.0, 0.0, 0.0), geom, goal_changed, 1.0)
    @test decision.action == :replan
    @test decision.reason == :goal_changed

    same_goal = SM.GuidanceHooks.RPOReplanningConfig(
        desired_goal_rtn=SVector{3, Float64}(2.0, 0.0, 0.0),
        goal_change_tolerance_m=0.05,
        safe_distance_m=0.1,
        retime_clearance_m=0.5,
    )
    decision = SM.GuidanceHooks.rpo_replanning_decision(plan, SVector{3, Float64}(0.0, 0.0, 0.0), geom, same_goal, 1.0)
    @test decision.action == :none
    @test decision.reason == :no_active_spheres
end

@testset "RPO planner comparison smoke csv-only outputs" begin
    points = [
        0.0 0.0;
        0.0 2.0;
        0.0 0.0
    ]
    geom = SM.RPOReferenceGeometry(
        SM.RPOStationGeometry(points; keepout_radius_m=0.2);
        chaser=SM.RPOCubeSatGeometry(dims_m=(0.1, 0.1, 0.1)),
    )
    pso_cfg = SM.rpo_740_mpc_final_pso_config(
        safe_distance_m=0.2;
        n_particles=8,
        n_iters=2,
        n_waypoints=1,
        adaptive_enable=false,
        refinement_enable=false,
        reexplore_enable=false,
        cull_enable=false,
        schedule_enable=false,
        sample_ds_m=0.5,
        retime_dt_s=0.5,
        retime_a_max_mps2=0.05,
        retime_max_steps=5_000,
    )
    @test pso_cfg.mass_kg == 5.0
    @test pso_cfg.retime_reaction_time_s == 0.25
    @test pso_cfg.retime_speed_scale == 0.5
    @test pso_cfg.retime_max_speed_mps == 0.25
    @test pso_cfg.search_margin_m == 10.0
    @test pso_cfg.reexplore_search_margin_scale == 2.0

    tracking = SM.RPOLQMPCTrackingSettings(
        dt_s=0.5,
        horizon=4,
        propellant_mass_kg=0.5,
        settle_time_s=2.0,
        final_position_tol_m=0.75,
        u_max_mps2=SVector{3, Float64}(0.05, 0.05, 0.05),
    )
    cfg = SM.RPOPlannerComparisonConfig(
        pso_config=pso_cfg,
        rrt_connect=SM.RPORRTConnectSettings(n_iters=25, shortcut_iters=4),
        rrt_star=SM.RPORRTStarSettings(n_iters=25, shortcut_iters=4),
        chomp=SM.RPOCHOMPSettings(n_iters=2),
        stomp=SM.RPOSTOMPSettings(n_iters=2, n_rollouts=3),
        optimizer=SM.RPOTrajectoryOptimizerSettings(match_hypr_runtime=false, force_full_iters=true),
        tracking=tracking,
        safe_distance_m=0.2,
        rng_seed=740,
    )
    cases = [
        SM.RPOPlannerComparisonCase(
            start_rtn=SVector{3, Float64}(-4.0, 1.0, 0.0),
            goal_rtn=SVector{3, Float64}(4.0, 1.0, 0.0),
            label="clear-offset",
        ),
    ]

    mktempdir() do output_dir
        cfg = SM.RPOPlannerComparisonConfig(
            pso_config=cfg.pso_config,
            rrt_connect=cfg.rrt_connect,
            rrt_star=cfg.rrt_star,
            chomp=cfg.chomp,
            stomp=cfg.stomp,
            optimizer=cfg.optimizer,
            tracking=cfg.tracking,
            safe_distance_m=cfg.safe_distance_m,
            output_dir=output_dir,
            write_plotly_outputs=false,
            rng_seed=cfg.rng_seed,
        )
        batch = SM.rpo_run_planner_comparison_batch(cases, geom, cfg)
        rows = SM.GuidanceHooks.rpo_flatten_planner_results(batch)
        @test SM.normalize_rpo_comparison_planner_type(:rrtconnect) == :rrt_connect
        @test SM.normalize_rpo_comparison_planner_type(:rrt_connect_bezier) == :rrt_connect_bezier
        @test SM.normalize_rpo_comparison_planner_type(:unrefined_pso) == :pso_unrefined
        @test SM.normalize_rpo_comparison_planner_type(:rrtstar) == :rrt_star
        @test Set([row.planner for row in rows]) == Set([:hypr, :pso_unrefined, :rrt_connect, :rrt_connect_bezier, :rrt_star, :chomp, :stomp])
        @test all(isfinite(row.fuel_used) && row.fuel_used > 0.0 for row in rows)
        @test all(isfinite(row.fuel_used_pct) && row.fuel_used_pct > 0.0 for row in rows)
        @test all(row.fuel_used_pct ≈ 100.0 * row.fuel_used / tracking.propellant_mass_kg for row in rows)
        @test all(isfinite(row.planner_compute_time) for row in rows)
        @test all(row.success for row in rows)
        @test SM.GuidanceHooks.rpo_comparison_metric_includes_result((success=false,), :success_pct)
        @test SM.GuidanceHooks.rpo_comparison_metric_includes_result((success=false,), :keepout_violations)
        @test !SM.GuidanceHooks.rpo_comparison_metric_includes_result((success=false,), :fuel_used_pct)

        metric_plot = SM.GuidanceHooks.rpo_comparison_metric_summary_plot(batch)
        box_traces = [trace for trace in metric_plot.data if get(trace.fields, :type, "") == "box"]
        @test !isempty(box_traces)
        @test all(get(trace.fields, :boxmean, false) === true for trace in box_traces)

        path_plot = SM.GuidanceHooks.rpo_comparison_path_family_plot(batch; planner=:hypr)
        trace_names = [String(get(trace.fields, :name, "")) for trace in path_plot.data]
        @test length(path_plot.data) == 5
        @test any(endswith(name, " desired") for name in trace_names)
        @test any(endswith(name, " tracked") for name in trace_names)
        @test !any(endswith(name, " path") for name in trace_names)
        @test !any(endswith(name, " retimed") for name in trace_names)

        cost_plot = SM.GuidanceHooks.rpo_comparison_cost_iteration_plot(batch; planner=:hypr)
        @test length(cost_plot.data) == 1

        outputs = SM.rpo_write_planner_comparison_outputs(batch)
        @test isfile(outputs.csv)
        @test outputs.plotly_outputs === false
        @test isnothing(outputs.metrics_plot)
        @test isnothing(outputs.path_plot)
        @test isempty(outputs.method_path_plots)
        @test isempty(outputs.cost_iteration_plots)
        @test filesize(outputs.csv) > 0
    end
end

@testset "RPO six-axis thrust allocation" begin
    thrusters = SM.SixAxisThrusterModel(max_thrust_n=SVector{6, Float64}(1, 1, 1, 1, 1, 1))
    forces = SM.rpo_allocate_six_axis_thrusters(SVector{3, Float64}(0.5, -0.25, 0.1), thrusters)
    @test forces ≈ SVector{6, Float64}(0.5, 0.0, 0.0, 0.25, 0.1, 0.0)
    body_force, body_torque = SM.rpo_thruster_wrench_body(forces, thrusters)
    @test body_force ≈ SVector{3, Float64}(0.5, -0.25, 0.1)
    @test body_torque ≈ SVector{3, Float64}(0.0, 0.0, 0.0)
end
