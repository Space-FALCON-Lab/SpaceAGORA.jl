"""Settings for LQ-MPC path tracking metrics in RPO planner comparisons."""
Base.@kwdef struct RPOLQMPCTrackingSettings
    dt_s::Float64 = 0.1
    mean_motion_radps::Float64 = 0.0011
    horizon::Int = 60
    mass_kg::Float64 = 5.0
    propellant_mass_kg::Float64 = 0.2
    isp_s::Float64 = 60.0
    g0_mps2::Float64 = 9.80665
    u_max_mps2::SVector{3, Float64} = SVector{3, Float64}(0.0125, 0.0125, 0.0125)
    q_pos::Float64 = 10.0
    q_vel::Float64 = 1.0
    r_accel::Float64 = 0.1
    qf_pos::Float64 = 50.0
    qf_vel::Float64 = 5.0
    settle_time_s::Float64 = 20.0
    final_position_tol_m::Float64 = 0.25
end

"""Single start-goal case definition for planner comparison batches."""
Base.@kwdef struct RPOPlannerComparisonCase
    start_rtn::SVector{3, Float64}
    goal_rtn::SVector{3, Float64}
    label::String = "case"
end

"""Configuration for running and exporting RPO planner comparison batches."""
const RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M = 0.5

Base.@kwdef struct RPOPlannerComparisonConfig
    planners::Vector{Symbol} = [:hypr, :pso_unrefined, :rrt_connect, :rrt_connect_bezier, :rrt_star, :chomp, :stomp]
    pso_config::RPOPSOConfig = rpo_740_mpc_final_pso_config(safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M)
    rrt_connect::RPORRTConnectSettings = RPORRTConnectSettings()
    rrt_star::RPORRTStarSettings = RPORRTStarSettings()
    chomp::RPOCHOMPSettings = RPOCHOMPSettings(n_iters=pso_config.n_iters)
    stomp::RPOSTOMPSettings = RPOSTOMPSettings(n_iters=pso_config.n_iters)
    optimizer::RPOTrajectoryOptimizerSettings = RPOTrajectoryOptimizerSettings()
    tracking::RPOLQMPCTrackingSettings = RPOLQMPCTrackingSettings()
    safe_distance_m::Float64 = pso_config.safe_distance_m > 0.0 ? pso_config.safe_distance_m : RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M
    output_dir::String = joinpath(pwd(), "rpo_planner_comparison")
    write_plotly_outputs::Bool = true
    write_failed_path_outputs::Bool = true
    rng_seed::Int = 740
    show_progress::Bool = true
end

"""Return the comparison config with planner cost safety pinned to the published 0.5 m value."""
function _rpo_comparison_config_with_fixed_safe_distance(cfg::RPOPlannerComparisonConfig)
    return RPOPlannerComparisonConfig(
        planners=cfg.planners,
        pso_config=rpo_pso_config(cfg.pso_config; safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M),
        rrt_connect=cfg.rrt_connect,
        rrt_star=cfg.rrt_star,
        chomp=cfg.chomp,
        stomp=cfg.stomp,
        optimizer=cfg.optimizer,
        tracking=cfg.tracking,
        safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M,
        output_dir=cfg.output_dir,
        write_plotly_outputs=cfg.write_plotly_outputs,
        write_failed_path_outputs=cfg.write_failed_path_outputs,
        rng_seed=cfg.rng_seed,
        show_progress=cfg.show_progress,
    )
end

"""Render a fixed-width text progress bar for planner comparison output."""
function _rpo_comparison_progress_bar(frac::Real; width::Integer=30)
    f = clamp(Float64(frac), 0.0, 1.0)
    filled = Int(floor(f * width))
    return "[" * repeat("=", filled) * repeat(" ", width - filled) * "]"
end

"""Print or update one planner comparison progress line."""
function _rpo_comparison_progress_line!(
    completed::Integer,
    total::Integer;
    planner::Symbol=:unknown,
    case_label::AbstractString="",
)
    frac = total <= 0 ? 1.0 : completed / total
    percent = lpad(string(round(100.0 * clamp(frac, 0.0, 1.0); digits=1)), 5)
    detail = planner == :unknown ? "" : "  $(rpo_comparison_planner_label(planner))"
    isempty(case_label) || (detail *= " / $(case_label)")
    line = "RPO planner comparison $(_rpo_comparison_progress_bar(frac)) $(percent)%  $(completed)/$(total)$(detail)"
    print("\r", rpad(line, 120))
    flush(stdout)
    completed >= total && println()
    return nothing
end

"""Build the tuned RPO HYPR configuration used by the 740 MPC comparison case."""
function rpo_740_mpc_final_pso_config(; safe_distance_m::Real=0.0, goal_collision_margin_m::Real=0.0, kwargs...)
    cfg = RPOPSOConfig(
        n_waypoints=5,
        n_particles=100,
        n_iters=60,
        w_len=1.0,
        w_obs=1.0e6,
        w_fuel=1.0,
        w_inertia=0.7,
        c1=1.4,
        c2=1.4,
        spread_scale=0.2,
        search_margin_m=10.0,
        sample_ds_m=0.05,
        curve_type=:bezier,
        cost_ref_distance_m=20.0,
        mass_kg=5.0,
        tf_s=120.0,
        isp_s=60.0,
        g0_mps2=9.80665,
        retime_dt_s=0.1,
        retime_reaction_time_s=0.25,
        retime_a_max_mps2=0.00625,
        retime_speed_scale=0.5,
        retime_min_speed_mps=0.0,
        retime_max_speed_mps=0.25,
        retime_max_steps=20_000,
        safe_distance_m=Float64(safe_distance_m),
        goal_collision_margin_m=Float64(goal_collision_margin_m),
        adaptive_enable=true,
        adaptive_allow_downscale=true,
        adaptive_complexity_weight=0.35,
        adaptive_n_waypoints_min=3,
        adaptive_n_waypoints_max=8,
        adaptive_n_particles_min=60,
        adaptive_n_particles_max=160,
        adaptive_n_iters_min=10,
        adaptive_n_iters_max=60,
        adaptive_w_len_min=0.5,
        adaptive_w_len_max=1.5,
        adaptive_w_obs_min=1.0e6,
        adaptive_w_obs_max=1.0e6,
        adaptive_w_inertia_min=0.4,
        adaptive_w_inertia_max=0.75,
        adaptive_c1_min=1.2,
        adaptive_c1_max=1.8,
        adaptive_c2_min=1.2,
        adaptive_c2_max=2.2,
        adaptive_spread_scale_min=0.05,
        adaptive_spread_scale_max=0.5,
        early_stopping_enable=true,
        early_stopping_patience=10,
        early_stopping_min_iters=15,
        early_stopping_min_rel_improvement=1.0e-4,
        early_stopping_require_feasible=true,
        cull_enable=true,
        cull_fraction_max=0.25,
        cull_start_iter=50,
        cull_noise_scale=0.3,
        cull_arc_velocity_scale=0.12,
        schedule_enable=true,
        schedule_w_end_fraction=0.65,
        schedule_c1_end_fraction=0.75,
        schedule_c2_end_fraction=1.25,
        schedule_w_min=0.25,
        schedule_c_min=0.5,
        schedule_c_max=2.5,
        stagnation_learning_enable=true,
        stagnation_learning_threshold=8,
        stagnation_learning_elite_fraction=0.10,
        stagnation_learning_max_blocks=2,
        probe_enable=true,
        probe_max_depth=2,
        probe_candidates=24,
        probe_offset_scale=1.0,
        probe_sample_ds_m=0.25,
        probe_seed=1234,
        reexplore_enable=true,
        reexplore_trigger_iter=10,
        reexplore_search_margin_scale=2.0,
        reexplore_waypoint_scale=1.5,
        reexplore_waypoint_increment=2,
        reexplore_max_waypoints=max(8 + 4, Int(ceil(1.5 * 8))),
        rrt_warmstart_enable=true,
        refinement_enable=true,
        refinement_start_iter=60,
        refinement_period=10,
        refinement_sample_ds_m=0.025,
        refinement_merge_distance_m=1.0,
        refinement_waypoint_passes=24,
        refinement_rounds=16,
        refinement_min_abs_cost_improvement=1.0e-8,
        refinement_min_rel_cost_improvement=1.0e-5,
        refinement_insert_straight_waypoints=true,
        refinement_straight_max_segment_length_m=2.0,
        refinement_straight_max_inserted=12,
        refinement_straight_clearance_margin_m=0.0,
    )
    return rpo_pso_config(cfg; kwargs...)
end

"""Build RRT-Connect settings for comparison planners with a caller-selected iteration cap."""
function _rpo_comparison_rrt_connect_settings(cfg::RPOPlannerComparisonConfig, n_iters::Integer)
    src = cfg.rrt_connect
    return RPORRTConnectSettings(
        n_iters=Int(n_iters),
        step_size_m=src.step_size_m,
        goal_sample_rate=src.goal_sample_rate,
        collision_sample_ds_m=src.collision_sample_ds_m,
        adaptive_collision_sampling_enable=src.adaptive_collision_sampling_enable,
        collision_max_sample_ds_m=src.collision_max_sample_ds_m,
        collision_far_clearance_m=src.collision_far_clearance_m,
        collision_sampling_power=src.collision_sampling_power,
        collision_safe_distance_fraction=src.collision_safe_distance_fraction,
        collision_obstacle_guard_fraction=src.collision_obstacle_guard_fraction,
        connect_max_steps=src.connect_max_steps,
        shortcut_iters=src.shortcut_iters,
    )
end

"""Build RRT* settings for comparison planners with a caller-selected iteration cap."""
function _rpo_comparison_rrt_star_settings(cfg::RPOPlannerComparisonConfig, n_iters::Integer)
    src = cfg.rrt_star
    return RPORRTStarSettings(
        n_iters=Int(n_iters),
        step_size_m=src.step_size_m,
        goal_sample_rate=src.goal_sample_rate,
        collision_sample_ds_m=src.collision_sample_ds_m,
        adaptive_collision_sampling_enable=src.adaptive_collision_sampling_enable,
        collision_max_sample_ds_m=src.collision_max_sample_ds_m,
        collision_far_clearance_m=src.collision_far_clearance_m,
        collision_sampling_power=src.collision_sampling_power,
        collision_safe_distance_fraction=src.collision_safe_distance_fraction,
        collision_obstacle_guard_fraction=src.collision_obstacle_guard_fraction,
        neighbor_radius_m=src.neighbor_radius_m,
        shortcut_iters=src.shortcut_iters,
    )
end

"""Generate the RRT-Connect seed used to warm-start trajectory optimizers."""
function _rpo_comparison_rrt_connect_seed_path(
    start_rtn,
    goal_rtn,
    geometry,
    base_cfg::RPOPSOConfig,
    cfg::RPOPlannerComparisonConfig,
    n_iters::Integer;
    runtime_limit_s::Real,
    rng,
)
    settings = _rpo_comparison_rrt_connect_settings(cfg, n_iters)
    plan = rpo_rrt_connect_plan_path(
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=cfg.safe_distance_m,
        settings=settings,
        max_runtime_s=runtime_limit_s,
        rng=rng,
    )
    return plan.path
end

"""Run one planner on one RPO comparison case and collect comparable metrics."""
function rpo_plan_comparison_path(
    planner_type,
    start_rtn,
    goal_rtn,
    geometry,
    cfg::RPOPlannerComparisonConfig;
    rng=Random.default_rng(),
    runtime_limit_s::Real=Inf,
)
    cfg = _rpo_comparison_config_with_fixed_safe_distance(cfg)
    planner = normalize_rpo_comparison_planner_type(planner_type)
    base_cfg = rpo_pso_config(cfg.pso_config; safe_distance_m=cfg.safe_distance_m)
    # Choose optimizer iterations from the runtime-limited override, HYPR matching, or planner default.
    runtime_limited_iters(default_iters) = isfinite(Float64(runtime_limit_s)) && cfg.optimizer.runtime_max_iters > 0 ?
        cfg.optimizer.runtime_max_iters :
        (cfg.optimizer.match_hypr_iters ? base_cfg.n_iters : default_iters)
    t0 = time_ns()
    if planner == :hypr
        plan = rpo_pso_plan_path(start_rtn, goal_rtn, geometry, base_cfg; safe_distance_m=cfg.safe_distance_m, rng=rng)
        elapsed = (time_ns() - t0) / 1.0e9
        iterations = length(plan.cost_history)
        return merge(plan, (
            planner_type=:hypr,
            planner_label=rpo_comparison_planner_label(:hypr),
            planner_compute_time=elapsed,
            planner_plan_time=elapsed,
            planner_refinement_time=0.0,
            planner_iteration_count=iterations,
        ))
    elseif planner == :pso_unrefined
        unrefined_cfg = rpo_pso_config(base_cfg; refinement_enable=false)
        plan = rpo_pso_plan_path(start_rtn, goal_rtn, geometry, unrefined_cfg; safe_distance_m=cfg.safe_distance_m, rng=rng)
        elapsed = (time_ns() - t0) / 1.0e9
        iterations = length(plan.cost_history)
        return merge(plan, (
            planner_type=:pso_unrefined,
            planner_label=rpo_comparison_planner_label(:pso_unrefined),
            planner_compute_time=elapsed,
            planner_plan_time=elapsed,
            planner_refinement_time=0.0,
            planner_iteration_count=iterations,
        ))
    elseif planner == :rrt_connect
        settings = _rpo_comparison_rrt_connect_settings(cfg, runtime_limited_iters(cfg.rrt_connect.n_iters))
        plan = rpo_rrt_connect_plan_path(
            start_rtn,
            goal_rtn,
            geometry,
            base_cfg;
            safe_distance_m=cfg.safe_distance_m,
            settings=settings,
            max_runtime_s=runtime_limit_s,
            rng=rng,
        )
        elapsed = (time_ns() - t0) / 1.0e9
        return merge(plan, (
            planner_type=:rrt_connect,
            planner_label=rpo_comparison_planner_label(:rrt_connect),
            planner_compute_time=elapsed,
            planner_plan_time=elapsed,
            planner_refinement_time=0.0,
            planner_iteration_count=plan.iterations,
        ))
    elseif planner == :rrt_connect_bezier
        settings = _rpo_comparison_rrt_connect_settings(cfg, runtime_limited_iters(cfg.rrt_connect.n_iters))
        plan = rpo_rrt_connect_bezier_plan_path(
            start_rtn,
            goal_rtn,
            geometry,
            base_cfg;
            safe_distance_m=cfg.safe_distance_m,
            settings=settings,
            max_runtime_s=runtime_limit_s,
            rng=rng,
        )
        elapsed = (time_ns() - t0) / 1.0e9
        return merge(plan, (
            planner_type=:rrt_connect_bezier,
            planner_label=rpo_comparison_planner_label(:rrt_connect_bezier),
            planner_compute_time=elapsed,
            planner_plan_time=elapsed,
            planner_refinement_time=0.0,
            planner_iteration_count=plan.iterations,
        ))
    elseif planner == :rrt_star
        settings = _rpo_comparison_rrt_star_settings(cfg, runtime_limited_iters(cfg.rrt_star.n_iters))
        plan = rpo_rrt_star_plan_path(
            start_rtn,
            goal_rtn,
            geometry,
            base_cfg;
            safe_distance_m=cfg.safe_distance_m,
            settings=settings,
            max_runtime_s=runtime_limit_s,
            rng=rng,
        )
        elapsed = (time_ns() - t0) / 1.0e9
        return merge(plan, (
            planner_type=:rrt_star,
            planner_label=rpo_comparison_planner_label(:rrt_star),
            planner_compute_time=elapsed,
            planner_plan_time=elapsed,
            planner_refinement_time=0.0,
            planner_iteration_count=plan.iterations,
        ))
    elseif planner == :chomp
        settings = RPOCHOMPSettings(
            n_iters=runtime_limited_iters(cfg.chomp.n_iters),
            learning_rate=cfg.chomp.learning_rate,
            gradient_eps=cfg.chomp.gradient_eps,
            w_smooth=cfg.chomp.w_smooth,
        )
        initial_path = _rpo_comparison_rrt_connect_seed_path(
            start_rtn,
            goal_rtn,
            geometry,
            base_cfg,
            cfg,
            runtime_limited_iters(cfg.rrt_connect.n_iters);
            runtime_limit_s=runtime_limit_s,
            rng=rng,
        )
        plan = rpo_chomp_plan_path(
            start_rtn,
            goal_rtn,
            geometry,
            base_cfg;
            safe_distance_m=cfg.safe_distance_m,
            settings=settings,
            optimizer=cfg.optimizer,
            max_runtime_s=runtime_limit_s,
            initial_path=initial_path,
        )
        elapsed = (time_ns() - t0) / 1.0e9
        return merge(plan, (
            planner_type=:chomp,
            planner_label=rpo_comparison_planner_label(:chomp),
            planner_compute_time=elapsed,
            planner_plan_time=elapsed,
            planner_refinement_time=0.0,
            planner_iteration_count=plan.iterations,
        ))
    else
        settings = RPOSTOMPSettings(
            n_iters=runtime_limited_iters(cfg.stomp.n_iters),
            n_rollouts=cfg.stomp.n_rollouts,
            noise_std=cfg.stomp.noise_std,
            lambda=cfg.stomp.lambda,
            update_step=cfg.stomp.update_step,
            w_smooth=cfg.stomp.w_smooth,
        )
        initial_path = _rpo_comparison_rrt_connect_seed_path(
            start_rtn,
            goal_rtn,
            geometry,
            base_cfg,
            cfg,
            runtime_limited_iters(cfg.rrt_connect.n_iters);
            runtime_limit_s=runtime_limit_s,
            rng=rng,
        )
        plan = rpo_stomp_plan_path(
            start_rtn,
            goal_rtn,
            geometry,
            base_cfg;
            safe_distance_m=cfg.safe_distance_m,
            settings=settings,
            optimizer=cfg.optimizer,
            max_runtime_s=runtime_limit_s,
            rng=rng,
            initial_path=initial_path,
        )
        elapsed = (time_ns() - t0) / 1.0e9
        return merge(plan, (
            planner_type=:stomp,
            planner_label=rpo_comparison_planner_label(:stomp),
            planner_compute_time=elapsed,
            planner_plan_time=elapsed,
            planner_refinement_time=0.0,
            planner_iteration_count=plan.iterations,
        ))
    end
end

"""Build a finite-horizon RPO reference preview for LQ-MPC tracking."""
function rpo_lqmpc_reference_preview(r_ref, v_ref, start_idx::Integer, horizon::Integer)
    out = zeros(6, horizon + 1)
    n_ref = size(r_ref, 2)
    for j in 0:horizon
        idx = min(start_idx + j, n_ref)
        out[1:3, j + 1] .= r_ref[:, idx]
        out[4:6, j + 1] .= v_ref[:, idx]
    end
    return out
end

"""Convert consumed propellant mass into a percent of configured initial mass."""
function rpo_lqmpc_tracking_fuel_used_pct(fuel_used_kg::Real, tracking::RPOLQMPCTrackingSettings)::Float64
    propellant_mass = Float64(tracking.propellant_mass_kg)
    return isfinite(propellant_mass) && propellant_mass > 0.0 ?
        100.0 * Float64(fuel_used_kg) / propellant_mass :
        NaN
end

"""Track a retimed RPO path with LQ-MPC and report tracking/fuel metrics."""
function rpo_track_retimed_path_lqmpc(path_rtn, goal_rtn, geometry, pso_cfg::RPOPSOConfig, tracking::RPOLQMPCTrackingSettings; safe_distance_m::Real=0.0)
    retime_cfg = rpo_pso_config(
        pso_cfg;
        retime_dt_s=tracking.dt_s,
        mass_kg=tracking.mass_kg,
        isp_s=tracking.isp_s,
        g0_mps2=tracking.g0_mps2,
    )
    t_ref, r_ref, v_ref = rpo_reference_from_path(path_rtn, geometry, retime_cfg; safe_distance_m=safe_distance_m)
    Q = Diagonal(vcat(fill(tracking.q_pos, 3), fill(tracking.q_vel, 3)))
    R = Diagonal(fill(tracking.r_accel, 3))
    Qf = Diagonal(vcat(fill(tracking.qf_pos, 3), fill(tracking.qf_vel, 3)))
    u_max = Vector{Float64}(tracking.u_max_mps2)
    u_min = -u_max
    ctrl = _control_module().init_rpo_lqmpc(
        tracking.mean_motion_radps,
        tracking.dt_s,
        Q,
        R,
        Qf,
        tracking.horizon;
        u_min=u_min,
        u_max=u_max,
    )
    x = zeros(6)
    x[1:3] .= r_ref[:, 1]
    n_plan_steps = max(size(r_ref, 2) - 1, 0)
    settle_steps = max(0, Int(ceil(tracking.settle_time_s / tracking.dt_s)))
    total_steps = n_plan_steps + settle_steps
    x_hist = zeros(6, total_steps + 1)
    u_hist = zeros(3, total_steps)
    x_hist[:, 1] .= x

    fuel_used = 0.0
    control_effort = 0.0
    saturated = 0
    min_clearance = Inf
    keepout_violations = 0

    for k in 1:total_steps
        ref_idx = min(k, size(r_ref, 2))
        x_ref = rpo_lqmpc_reference_preview(r_ref, v_ref, ref_idx, tracking.horizon)
        u = _control_module().rpo_lqmpc_control(ctrl, x, x_ref)
        u_vec = Vector{Float64}(u)
        u_norm = norm(u_vec)
        fuel_used += tracking.mass_kg * sum(abs, u_vec) * tracking.dt_s / max(tracking.isp_s * tracking.g0_mps2, 1.0e-9)
        control_effort += u_norm * tracking.dt_s
        saturated += any(abs.(u_vec) .>= (0.999 .* u_max)) ? 1 : 0
        x .= ctrl.Ad * x .+ ctrl.Bd * u_vec
        x_hist[:, k + 1] .= x
        u_hist[:, k] .= u_vec

        clearance = rpo_clearance_distance_to_station(x[1:3], geometry)
        min_clearance = min(min_clearance, clearance)
        keepout_violations += clearance < 0.0 ? 1 : 0
    end

    final_error = norm(x[1:3] - Vector{Float64}(goal_rtn))
    success = final_error <= tracking.final_position_tol_m && min_clearance + 1.0e-9 >= 0.0
    return (
        success=success,
        fuel_used=fuel_used,
        fuel_used_pct=rpo_lqmpc_tracking_fuel_used_pct(fuel_used, tracking),
        control_effort_total=control_effort,
        translational_control_effort=control_effort,
        thrust_saturation_fraction=total_steps > 0 ? saturated / total_steps : 0.0,
        min_clearance=min_clearance,
        keepout_violations=keepout_violations,
        final_pos_error=final_error,
        planned_travel_duration=isempty(t_ref) ? 0.0 : t_ref[end],
        actual_travel_duration=tracking.dt_s * total_steps,
        t_ref=t_ref,
        r_ref_rtn=r_ref,
        v_ref_rtn=v_ref,
        x_hist=x_hist,
        u_hist=u_hist,
    )
end

"""Run all configured RPO planners across all comparison cases."""
function rpo_run_planner_comparison_batch(cases, geometry, cfg::RPOPlannerComparisonConfig=RPOPlannerComparisonConfig())
    cfg = _rpo_comparison_config_with_fixed_safe_distance(cfg)
    comparison_cases = collect(cases)
    planner_types = [normalize_rpo_comparison_planner_type(p) for p in cfg.planners]
    if cfg.optimizer.match_hypr_runtime && (:hypr in planner_types)
        planner_types = vcat([:hypr], [p for p in planner_types if p != :hypr])
    end

    results_by_planner = Dict{Symbol, Vector{NamedTuple}}()
    plans_by_planner = Dict{Symbol, Vector{NamedTuple}}()
    master_rng = MersenneTwister(cfg.rng_seed)
    total_runs = length(planner_types) * length(comparison_cases)
    completed_runs = 0
    cfg.show_progress && _rpo_comparison_progress_line!(completed_runs, total_runs)

    for planner in planner_types
        planner_idx = findfirst(==(planner), planner_types)
        planner_results = NamedTuple[]
        planner_plans = NamedTuple[]
        for (case_idx, case) in enumerate(comparison_cases)
            runtime_limit = cfg.optimizer.runtime_limit_s
            if planner != :hypr && cfg.optimizer.match_hypr_runtime && haskey(plans_by_planner, :hypr)
                runtime_limit = plans_by_planner[:hypr][case_idx].planner_plan_time
            end
            cfg.show_progress && _rpo_comparison_progress_line!(
                completed_runs,
                total_runs;
                planner=planner,
                case_label=case.label,
            )
            rng = MersenneTwister(rand(master_rng, UInt) + UInt(1000 * planner_idx + case_idx))
            plan = rpo_plan_comparison_path(
                planner,
                case.start_rtn,
                case.goal_rtn,
                geometry,
                cfg;
                rng=rng,
                runtime_limit_s=runtime_limit,
            )
            tracking = rpo_track_retimed_path_lqmpc(
                plan.path,
                case.goal_rtn,
                geometry,
                plan.config,
                cfg.tracking;
                safe_distance_m=cfg.safe_distance_m,
            )
            result = merge((
                    planner=planner,
                    planner_label=rpo_comparison_planner_label(planner),
                    case_id=case_idx,
                    case_label=case.label,
                    start_x=case.start_rtn[1],
                    start_y=case.start_rtn[2],
                    start_z=case.start_rtn[3],
                    goal_x=case.goal_rtn[1],
                    goal_y=case.goal_rtn[2],
                    goal_z=case.goal_rtn[3],
                    planner_compute_time=plan.planner_compute_time,
                    planner_plan_time=plan.planner_plan_time,
                    planner_refinement_time=plan.planner_refinement_time,
                    planner_iteration_count=plan.planner_iteration_count,
                    planner_cost=plan.cost,
                    refinement_improved=hasproperty(plan, :refinement_improved) ? plan.refinement_improved : false,
                ),
                (
                    success=tracking.success,
                    fuel_used=tracking.fuel_used,
                    fuel_used_pct=tracking.fuel_used_pct,
                    control_effort_total=tracking.control_effort_total,
                    translational_control_effort=tracking.translational_control_effort,
                    thrust_saturation_fraction=tracking.thrust_saturation_fraction,
                    min_clearance=tracking.min_clearance,
                    keepout_violations=tracking.keepout_violations,
                    final_pos_error=tracking.final_pos_error,
                    planned_travel_duration=tracking.planned_travel_duration,
                    actual_travel_duration=tracking.actual_travel_duration,
                ),
            )
            push!(planner_results, result)
            push!(planner_plans, merge(plan, (tracking=tracking, case=case)))
            completed_runs += 1
            cfg.show_progress && _rpo_comparison_progress_line!(
                completed_runs,
                total_runs;
                planner=planner,
                case_label=case.label,
            )
        end
        results_by_planner[planner] = planner_results
        plans_by_planner[planner] = planner_plans
    end

    return (
        cases=comparison_cases,
        geometry=geometry,
        planner_types=planner_types,
        results_by_planner=results_by_planner,
        plans_by_planner=plans_by_planner,
        config=cfg,
    )
end

"""Flatten nested planner-comparison output into one row per planner result."""
function rpo_flatten_planner_results(batch)
    rows = NamedTuple[]
    for planner in batch.planner_types
        append!(rows, batch.results_by_planner[planner])
    end
    return rows
end

"""Write named-tuple rows to a CSV file with stable column ordering."""
function rpo_write_namedtuple_csv(path::AbstractString, rows)
    isempty(rows) && return path
    mkpath(dirname(path))
    names = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(names), ","))
        for row in rows
            vals = String[]
            for name in names
                value = getproperty(row, name)
                push!(vals, replace(string(value), "," => ";"))
            end
            println(io, join(vals, ","))
        end
    end
    return path
end

"""Compute mean metric values grouped by planner label."""
function rpo_group_metric_mean(results, metric::Symbol)
    values = Float64[]
    for row in results
        value = getproperty(row, metric)
        value isa Bool && (value = value ? 1.0 : 0.0)
        f = Float64(value)
        isfinite(f) && push!(values, f)
    end
    isempty(values) && return NaN
    return sum(values) / length(values)
end

"""Return plot/export metadata for RPO planner comparison metrics."""
function rpo_comparison_metric_specs()
    return [
        (:success_pct, "Success", "success (%)"),
        (:planner_compute_time, "Planner Runtime", "s"),
        (:fuel_used_pct, "Fuel", "%"),
        (:control_effort_total, "Control Effort", "m/s"),
        (:thrust_saturation_fraction, "Thruster Sat.", "fraction"),
        (:min_clearance, "Min Clearance", "m"),
        (:final_pos_error, "Goal Error", "m"),
        (:actual_travel_duration, "Travel Duration", "s"),
        (:keepout_violations, "Keepout Violations", "count"),
    ]
end

"""Read one named comparison metric from a flattened result row."""
function rpo_comparison_metric_value(row, metric::Symbol)
    metric == :success_pct && return row.success ? 100.0 : 0.0
    value = getproperty(row, metric)
    value isa Bool && return value ? 1.0 : 0.0
    return Float64(value)
end

"""Decide whether a planner result contributes to a given metric aggregate."""
@inline rpo_comparison_metric_includes_result(row, metric::Symbol) = metric in (:success_pct, :keepout_violations) || row.success

"""Return the Plotly axis name for a metric subplot index."""
function rpo_comparison_axis_name(prefix::AbstractString, idx::Integer)
    idx == 1 && return Symbol(prefix)
    return Symbol(prefix, idx)
end

"""Return the Plotly trace axis reference for a metric subplot index."""
function rpo_comparison_trace_axis_name(prefix::AbstractString, idx::Integer)
    idx == 1 && return prefix
    return string(prefix, idx)
end

"""Create summary bar plots for planner comparison metrics."""
function rpo_comparison_metric_summary_plot(batch)
    specs = rpo_comparison_metric_specs()
    traces = PlotlyJS.GenericTrace[]
    rows = 3
    cols = 3
    layout_kwargs = Dict{Symbol, Any}(
        :title => "RPO Planner Comparison Metrics",
        :height => 920,
        :width => 1200,
        :showlegend => false,
        :margin => PlotlyJS.attr(l=70, r=30, t=80, b=50),
        :plot_bgcolor => "rgb(250,250,250)",
    )
    xgap = 0.055
    ygap = 0.095
    xspan = (1.0 - xgap * (cols - 1)) / cols
    yspan = (1.0 - ygap * (rows - 1)) / rows

    for (idx, (metric, label, units)) in enumerate(specs)
        row_idx = div(idx - 1, cols) + 1
        col_idx = mod(idx - 1, cols) + 1
        x0 = (col_idx - 1) * (xspan + xgap)
        x1 = x0 + xspan
        y1 = 1.0 - (row_idx - 1) * (yspan + ygap)
        y0 = y1 - yspan
        xaxis = rpo_comparison_trace_axis_name("x", idx)
        yaxis = rpo_comparison_trace_axis_name("y", idx)

        xvals = String[]
        yvals = Float64[]
        case_labels = String[]
        for planner in batch.planner_types
            planner_label = rpo_comparison_planner_label(planner)
            for result in batch.results_by_planner[planner]
                rpo_comparison_metric_includes_result(result, metric) || continue
                value = rpo_comparison_metric_value(result, metric)
                push!(xvals, planner_label)
                push!(yvals, value)
                push!(case_labels, "case $(result.case_label)")
            end
        end

        push!(traces, PlotlyJS.box(
            x=xvals,
            y=yvals,
            name=label,
            boxpoints=false,
            boxmean=true,
            fillcolor="rgba(64,120,165,0.28)",
            line=PlotlyJS.attr(color="rgb(42,86,125)", width=2),
            marker=PlotlyJS.attr(color="rgb(42,86,125)"),
            xaxis=xaxis,
            yaxis=yaxis,
            hovertemplate="$(label)<br>%{x}: %{y:.5g} $(units)<extra></extra>",
        ))
        push!(traces, PlotlyJS.scatter(
            x=xvals,
            y=yvals,
            mode="markers",
            name="$(label) cases",
            marker=PlotlyJS.attr(size=8, color="rgb(205,92,55)", line=PlotlyJS.attr(color="white", width=1)),
            text=case_labels,
            xaxis=xaxis,
            yaxis=yaxis,
            hovertemplate="%{text}<br>%{x}<br>$(label): %{y:.5g} $(units)<extra></extra>",
        ))

        layout_kwargs[rpo_comparison_axis_name("xaxis", idx)] = PlotlyJS.attr(
            domain=[x0, x1],
            anchor=yaxis,
            title=label,
        )
        layout_kwargs[rpo_comparison_axis_name("yaxis", idx)] = PlotlyJS.attr(
            domain=[y0, y1],
            anchor=xaxis,
            title=units,
            zeroline=true,
            gridcolor="rgb(225,225,225)",
        )
    end

    return PlotlyJS.Plot(traces, PlotlyJS.Layout(; layout_kwargs...))
end

"""Build a Plotly mesh trace for the station keepout geometry."""
function rpo_comparison_station_mesh_trace(triangles)
    tris = Matrix{Float64}(triangles)
    ntri = size(tris, 2) ÷ 3
    vertices = zeros(Float64, 3, 3 * ntri)
    i = zeros(Int, ntri)
    j = zeros(Int, ntri)
    k = zeros(Int, ntri)
    for tri_idx in 1:ntri
        src = 3 * (tri_idx - 1)
        vertices[:, src + 1] .= tris[:, src + 1]
        vertices[:, src + 2] .= tris[:, src + 2]
        vertices[:, src + 3] .= tris[:, src + 3]
        i[tri_idx] = src
        j[tri_idx] = src + 1
        k[tri_idx] = src + 2
    end
    return PlotlyJS.mesh3d(
        x=vertices[1, :],
        y=vertices[2, :],
        z=vertices[3, :],
        i=i,
        j=j,
        k=k,
        color="rgb(150,160,170)",
        opacity=0.55,
        flatshading=true,
        name="Gateway mesh",
        hoverinfo="skip",
    )
end

"""Build a downsampled Plotly station point-cloud trace."""
function rpo_comparison_station_trace(batch; max_points::Integer=2_000)
    if hasproperty(batch, :station_triangles) && batch.station_triangles !== nothing
        return rpo_comparison_station_mesh_trace(batch.station_triangles)
    end
    pts = batch.geometry.station.points_body
    stride = max(1, Int(ceil(size(pts, 2) / max(1, max_points))))
    idxs = 1:stride:size(pts, 2)
    return PlotlyJS.scatter3d(
        x=pts[1, idxs],
        y=pts[2, idxs],
        z=pts[3, idxs],
        mode="markers",
        marker=PlotlyJS.attr(size=2, color="rgb(150,160,170)", opacity=0.30),
        name="station point cloud",
        hoverinfo="skip",
    )
end

"""Create a 3D Plotly view of planner paths for all comparison cases."""
function rpo_comparison_path_family_plot(batch; planner=nothing)
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    planners = planner === nothing ? batch.planner_types : [normalize_rpo_comparison_planner_type(planner)]
    for planner_type in planners
        label = rpo_comparison_planner_label(planner_type)
        for plan in batch.plans_by_planner[planner_type]
            case_label = plan.case.label
            ref = plan.tracking.r_ref_rtn
            actual = plan.tracking.x_hist[1:3, :]
            push!(traces, PlotlyJS.scatter3d(
                x=ref[1, :],
                y=ref[2, :],
                z=ref[3, :],
                mode="lines",
                line=PlotlyJS.attr(width=4, dash="dot"),
                name="$(label) case $(case_label) desired",
            ))
            push!(traces, PlotlyJS.scatter3d(
                x=actual[1, :],
                y=actual[2, :],
                z=actual[3, :],
                mode="lines",
                line=PlotlyJS.attr(width=4),
                name="$(label) case $(case_label) tracked",
            ))
            push!(traces, PlotlyJS.scatter3d(
                x=[plan.case.start_rtn[1]],
                y=[plan.case.start_rtn[2]],
                z=[plan.case.start_rtn[3]],
                mode="markers",
                marker=PlotlyJS.attr(size=6, color="rgb(45,150,90)", symbol="circle"),
                name="$(label) case $(case_label) start",
            ))
            push!(traces, PlotlyJS.scatter3d(
                x=[plan.case.goal_rtn[1]],
                y=[plan.case.goal_rtn[2]],
                z=[plan.case.goal_rtn[3]],
                mode="markers",
                marker=PlotlyJS.attr(size=7, color="rgb(190,60,55)", symbol="diamond"),
                name="$(label) case $(case_label) goal",
            ))
        end
    end
    title = planner === nothing ? "RPO Planner Comparison Path Family" :
        "$(rpo_comparison_planner_label(planner)) RPO Path Family"
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title=title,
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis_title="radial (m)",
                yaxis_title="along-track (m)",
                zaxis_title="cross-track (m)",
            ),
        ),
    )
end

"""Create a 3D Plotly view for one planner/case path."""
function rpo_comparison_single_path_plot(batch, plan; planner)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    label = rpo_comparison_planner_label(planner_type)
    case_label = plan.case.label
    ref = plan.tracking.r_ref_rtn
    actual = plan.tracking.x_hist[1:3, :]
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    push!(traces, PlotlyJS.scatter3d(
        x=ref[1, :],
        y=ref[2, :],
        z=ref[3, :],
        mode="lines",
        line=PlotlyJS.attr(width=5, dash="dot", color="rgb(45,95,170)"),
        name="desired",
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=actual[1, :],
        y=actual[2, :],
        z=actual[3, :],
        mode="lines",
        line=PlotlyJS.attr(width=5, color="rgb(205,92,55)"),
        name="tracked",
    ))
    if hasproperty(plan, :raw_path)
        raw = Matrix{Float64}(plan.raw_path)
        push!(traces, PlotlyJS.scatter3d(
            x=raw[1, :],
            y=raw[2, :],
            z=raw[3, :],
            mode="lines+markers",
            line=PlotlyJS.attr(width=2, dash="dash", color="rgb(110,110,110)"),
            marker=PlotlyJS.attr(size=3, color="rgb(110,110,110)"),
            name="raw planner path",
        ))
    end
    push!(traces, PlotlyJS.scatter3d(
        x=[plan.case.start_rtn[1]],
        y=[plan.case.start_rtn[2]],
        z=[plan.case.start_rtn[3]],
        mode="markers",
        marker=PlotlyJS.attr(size=7, color="rgb(45,150,90)", symbol="circle"),
        name="start",
    ))
    push!(traces, PlotlyJS.scatter3d(
        x=[plan.case.goal_rtn[1]],
        y=[plan.case.goal_rtn[2]],
        z=[plan.case.goal_rtn[3]],
        mode="markers",
        marker=PlotlyJS.attr(size=8, color="rgb(190,60,55)", symbol="diamond"),
        name="goal",
    ))
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="$(label) Failed Case $(case_label) Path",
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis_title="radial (m)",
                yaxis_title="along-track (m)",
                zaxis_title="cross-track (m)",
            ),
        ),
    )
end

"""Create a 3D Plotly view containing all failed paths for one planner."""
function rpo_comparison_failed_paths_plot(batch; planner, failed_indices)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    label = rpo_comparison_planner_label(planner_type)
    traces = PlotlyJS.GenericTrace[rpo_comparison_station_trace(batch)]
    plans = batch.plans_by_planner[planner_type]
    for idx in failed_indices
        plan = plans[idx]
        case_label = plan.case.label
        ref = plan.tracking.r_ref_rtn
        actual = plan.tracking.x_hist[1:3, :]
        push!(traces, PlotlyJS.scatter3d(
            x=ref[1, :],
            y=ref[2, :],
            z=ref[3, :],
            mode="lines",
            line=PlotlyJS.attr(width=4, dash="dot"),
            name="case $(case_label) desired",
        ))
        push!(traces, PlotlyJS.scatter3d(
            x=actual[1, :],
            y=actual[2, :],
            z=actual[3, :],
            mode="lines",
            line=PlotlyJS.attr(width=4),
            name="case $(case_label) tracked",
        ))
        if hasproperty(plan, :raw_path)
            raw = Matrix{Float64}(plan.raw_path)
            push!(traces, PlotlyJS.scatter3d(
                x=raw[1, :],
                y=raw[2, :],
                z=raw[3, :],
                mode="lines+markers",
                line=PlotlyJS.attr(width=2, dash="dash"),
                marker=PlotlyJS.attr(size=3),
                name="case $(case_label) raw planner path",
            ))
        end
        push!(traces, PlotlyJS.scatter3d(
            x=[plan.case.start_rtn[1]],
            y=[plan.case.start_rtn[2]],
            z=[plan.case.start_rtn[3]],
            mode="markers",
            marker=PlotlyJS.attr(size=6, color="rgb(45,150,90)", symbol="circle"),
            name="case $(case_label) start",
        ))
        push!(traces, PlotlyJS.scatter3d(
            x=[plan.case.goal_rtn[1]],
            y=[plan.case.goal_rtn[2]],
            z=[plan.case.goal_rtn[3]],
            mode="markers",
            marker=PlotlyJS.attr(size=7, color="rgb(190,60,55)", symbol="diamond"),
            name="case $(case_label) goal",
        ))
    end
    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="$(label) Failed Case Paths ($(length(failed_indices)))",
            scene=PlotlyJS.attr(
                aspectmode="data",
                xaxis_title="radial (m)",
                yaxis_title="along-track (m)",
                zaxis_title="cross-track (m)",
            ),
        ),
    )
end

"""Choose x-axis values for a planner cost-history trace."""
function rpo_comparison_cost_iteration_xvalues(plan, costs)
    n = length(costs)
    n == 0 && return Int[]
    iterations = hasproperty(plan, :planner_iteration_count) ? Int(plan.planner_iteration_count) :
        (hasproperty(plan, :iterations) ? Int(plan.iterations) : n)
    if n == 1
        return [max(iterations, 0)]
    end
    return collect(1:n)
end

"""Create cost-versus-iteration traces for planners that report histories."""
function rpo_comparison_cost_iteration_plot(batch; planner)
    planner_type = normalize_rpo_comparison_planner_type(planner)
    haskey(batch.plans_by_planner, planner_type) ||
        throw(ArgumentError("No planner results found for $(planner_type)."))

    label = rpo_comparison_planner_label(planner_type)
    traces = PlotlyJS.GenericTrace[]
    for plan in batch.plans_by_planner[planner_type]
        hasproperty(plan, :cost_history) || continue
        costs = [Float64(c) for c in plan.cost_history if isfinite(Float64(c))]
        isempty(costs) && continue
        iters = rpo_comparison_cost_iteration_xvalues(plan, costs)
        push!(traces, PlotlyJS.scatter(
            x=iters,
            y=costs,
            mode=length(costs) == 1 ? "markers" : "lines+markers",
            name="case $(plan.case.label)",
            line=PlotlyJS.attr(width=2),
            marker=PlotlyJS.attr(size=7),
            hovertemplate="case $(plan.case.label)<br>iteration %{x}<br>cost %{y:.6g}<extra></extra>",
        ))
    end

    isempty(traces) && push!(traces, PlotlyJS.scatter(
        x=Int[],
        y=Float64[],
        mode="markers",
        name="no finite cost history",
    ))

    return PlotlyJS.Plot(
        traces,
        PlotlyJS.Layout(
            title="$(label) Cost vs Iteration",
            height=620,
            width=920,
            hovermode="closest",
            margin=PlotlyJS.attr(l=75, r=35, t=75, b=65),
            plot_bgcolor="rgb(250,250,250)",
            xaxis=PlotlyJS.attr(
                title="Iteration",
                zeroline=true,
                gridcolor="rgb(225,225,225)",
            ),
            yaxis=PlotlyJS.attr(
                title="Planner Cost",
                rangemode="tozero",
                gridcolor="rgb(225,225,225)",
            ),
        ),
    )
end

"""Return a filesystem-safe token for generated planner-comparison artifact names."""
function rpo_comparison_artifact_slug(value)
    token = replace(lowercase(String(value)), r"[^a-z0-9_=-]+" => "_")
    token = strip(token, '_')
    return isempty(token) ? "case" : token
end

"""Write one path HTML per planner containing all failed planner/case results."""
function rpo_write_failed_path_outputs(batch, output_dir::AbstractString)
    failed_dir = joinpath(output_dir, "rpo_failed_path_plots")
    failed_path_plots = Dict{Symbol, String}()
    for planner in batch.planner_types
        results = batch.results_by_planner[planner]
        failed_indices = [idx for idx in eachindex(results) if !results[idx].success]
        isempty(failed_indices) && continue
        mkpath(failed_dir)
        path = joinpath(failed_dir, "rpo_failed_paths_$(planner).html")
        PlotlyJS.savefig(rpo_comparison_failed_paths_plot(batch; planner=planner, failed_indices=failed_indices), path)
        failed_path_plots[planner] = path
    end
    return failed_path_plots
end

"""Write planner comparison CSVs and plots to the requested output directory."""
function rpo_write_planner_comparison_outputs(
    batch;
    output_dir::AbstractString=batch.config.output_dir,
    write_plotly_outputs::Bool=batch.config.write_plotly_outputs,
    write_failed_path_outputs::Bool=batch.config.write_failed_path_outputs,
)
    mkpath(output_dir)
    rows = rpo_flatten_planner_results(batch)
    csv_path = rpo_write_namedtuple_csv(joinpath(output_dir, "rpo_planner_comparison_results.csv"), rows)

    metric_plot_path = nothing
    path_plot_path = nothing
    method_path_plots = Dict{Symbol, String}()
    cost_iteration_plots = Dict{Symbol, String}()
    failed_path_plots = Dict{Symbol, String}()

    if write_plotly_outputs
        metric_plot_path = joinpath(output_dir, "rpo_planner_comparison_metrics.html")
        PlotlyJS.savefig(rpo_comparison_metric_summary_plot(batch), metric_plot_path)

        path_plot_path = joinpath(output_dir, "rpo_planner_comparison_paths.html")
        PlotlyJS.savefig(rpo_comparison_path_family_plot(batch), path_plot_path)

        for planner in batch.planner_types
            method_path = joinpath(output_dir, "rpo_planner_paths_$(planner).html")
            PlotlyJS.savefig(rpo_comparison_path_family_plot(batch; planner=planner), method_path)
            method_path_plots[planner] = method_path
        end

        for planner in batch.planner_types
            cost_path = joinpath(output_dir, "rpo_planner_cost_vs_iteration_$(planner).html")
            PlotlyJS.savefig(rpo_comparison_cost_iteration_plot(batch; planner=planner), cost_path)
            cost_iteration_plots[planner] = cost_path
        end
    end
    if write_failed_path_outputs
        failed_path_plots = rpo_write_failed_path_outputs(batch, output_dir)
    end

    return (
        csv=csv_path,
        metrics_plot=metric_plot_path,
        path_plot=path_plot_path,
        method_path_plots=method_path_plots,
        cost_iteration_plots=cost_iteration_plots,
        failed_path_plots=failed_path_plots,
        plotly_outputs=write_plotly_outputs,
        failed_path_outputs=write_failed_path_outputs,
    )
end
