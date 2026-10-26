include("two_body_plant.jl")

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
    plant::Symbol = :hcw
    actuator_frame::Symbol = :inertial
    minimum_throttle_fraction::Float64 = 0.0
    max_thrust_n::Float64 = 0.05
end

"""Single start-goal case definition for planner comparison batches."""
Base.@kwdef struct RPOPlannerComparisonCase
    start_rtn::SVector{3, Float64}
    goal_rtn::SVector{3, Float64}
    label::String = "case"
end

"""Safety distance (m) of the published planner comparison; the harness pins every planner's cost to it."""
const RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M = 0.5

"""
    RPOPlannerComparisonConfig(; planners, pso_config, rrt_connect, rrt_star, chomp, stomp, optimizer, tracking, ...)

Configuration for running and exporting RPO planner comparison batches.

`planners` lists the planner objects to compare, in order: `HYPRComparisonPlanner`,
`UnrefinedPSOPlanner`, `RRTConnectPlanner`, `RRTConnectBezierPlanner`,
`RRTStarPlanner`, `CHOMPPlanner` and `STOMPPlanner`. The default is all seven,
built from this configuration's `pso_config`, `rrt_connect`, `rrt_star`, `chomp`,
`stomp` and `optimizer`. An object carries its own settings, so a configuration
rebuilt from these fields with other settings keeps the objects' settings.
"""
struct RPOPlannerComparisonConfig
    planners::Vector{AbstractRPOComparisonPlanner}
    pso_config::RPOPSOConfig
    rrt_connect::RPORRTConnectSettings
    rrt_star::RPORRTStarSettings
    chomp::RPOCHOMPSettings
    stomp::RPOSTOMPSettings
    optimizer::RPOTrajectoryOptimizerSettings
    tracking::RPOLQMPCTrackingSettings
    safe_distance_m::Float64
    output_dir::String
    write_plotly_outputs::Bool
    write_failed_path_outputs::Bool
    rng_seed::Int
    show_progress::Bool
end

function RPOPlannerComparisonConfig(;
    pso_config::RPOPSOConfig=rpo_740_mpc_final_pso_config(safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M),
    rrt_connect::RPORRTConnectSettings=RPORRTConnectSettings(),
    rrt_star::RPORRTStarSettings=RPORRTStarSettings(),
    chomp::RPOCHOMPSettings=RPOCHOMPSettings(n_iters=pso_config.n_iters),
    stomp::RPOSTOMPSettings=RPOSTOMPSettings(n_iters=pso_config.n_iters),
    optimizer::RPOTrajectoryOptimizerSettings=RPOTrajectoryOptimizerSettings(),
    tracking::RPOLQMPCTrackingSettings=RPOLQMPCTrackingSettings(),
    safe_distance_m=pso_config.safe_distance_m > 0.0 ? pso_config.safe_distance_m : RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M,
    output_dir=joinpath(pwd(), "rpo_planner_comparison"),
    write_plotly_outputs=true,
    write_failed_path_outputs=true,
    rng_seed=740,
    show_progress=true,
    planners=_rpo_comparison_planners(pso_config, rrt_connect, rrt_star, chomp, stomp, optimizer),
)
    return RPOPlannerComparisonConfig(planners, pso_config, rrt_connect, rrt_star, chomp, stomp, optimizer,
        tracking, safe_distance_m, output_dir, write_plotly_outputs, write_failed_path_outputs, rng_seed,
        show_progress)
end

const _RPO_COMPARISON_PLANNER_IDS = (:hypr, :pso_unrefined, :rrt_connect, :rrt_connect_bezier, :rrt_star, :chomp, :stomp)

# The object a canonical identifier selects, built from configuration settings fields.
function _rpo_comparison_planner(id::Symbol, pso_config, rrt_connect, rrt_star, chomp, stomp, optimizer)
    id === :hypr && return HYPRComparisonPlanner(pso_config)
    id === :pso_unrefined && return UnrefinedPSOPlanner(pso_config)
    id === :rrt_connect && return RRTConnectPlanner(rrt_connect; config=pso_config, optimizer=optimizer)
    id === :rrt_connect_bezier && return RRTConnectBezierPlanner(rrt_connect; config=pso_config, optimizer=optimizer)
    id === :rrt_star && return RRTStarPlanner(rrt_star; config=pso_config, optimizer=optimizer)
    id === :chomp && return CHOMPPlanner(chomp; config=pso_config, optimizer=optimizer, seed=rrt_connect)
    id === :stomp && return STOMPPlanner(stomp; config=pso_config, optimizer=optimizer, seed=rrt_connect)
    throw(ArgumentError("Unsupported RPO comparison planner $(id)."))
end
_rpo_comparison_planners(pso_config, rrt_connect, rrt_star, chomp, stomp, optimizer) =
    AbstractRPOComparisonPlanner[_rpo_comparison_planner(id, pso_config, rrt_connect, rrt_star, chomp, stomp, optimizer)
        for id in _RPO_COMPARISON_PLANNER_IDS]

"""Return the comparison config with planner cost safety pinned to the published 0.5 m value."""
function _rpo_comparison_config_with_fixed_safe_distance(cfg::RPOPlannerComparisonConfig)
    return _with_fields(cfg;
        pso_config=rpo_pso_config(cfg.pso_config; safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M),
        safe_distance_m=RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M,
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

"""
Run one planner on one RPO comparison case and collect comparable metrics. The
required `rng` is the only random stream the planner draws from.

`planner` is one of the seven comparison planner objects, which runs with its own
settings and this comparison's fixed safety distance.
"""
function rpo_plan_comparison_path(
    planner::AbstractRPOComparisonPlanner,
    start_rtn,
    goal_rtn,
    geometry,
    cfg::RPOPlannerComparisonConfig;
    rng::Random.AbstractRNG,
    runtime_limit_s::Real=Inf,
    warmstart_path=nothing,
)
    cfg = _rpo_comparison_config_with_fixed_safe_distance(cfg)
    HYPRSupport.require_hypr()
    id = normalize_rpo_comparison_planner_type(planner)
    base_cfg = rpo_pso_config(planner.config; safe_distance_m=cfg.safe_distance_m)
    t0 = time_ns()
    plan, iterations = _rpo_comparison_planner_path(
        planner,
        start_rtn,
        goal_rtn,
        geometry,
        base_cfg;
        safe_distance_m=cfg.safe_distance_m,
        rng=rng,
        runtime_limit_s=runtime_limit_s,
        warmstart_path=warmstart_path,
    )
    elapsed = (time_ns() - t0) / 1.0e9
    return merge(plan, (
        planner_type=id,
        planner_label=rpo_comparison_planner_label(id),
        planner_compute_time=elapsed,
        planner_plan_time=elapsed,
        planner_refinement_time=0.0,
        planner_iteration_count=iterations,
    ))
end

"""Build a finite-horizon RPO reference preview for LQ-MPC tracking."""
function rpo_lqmpc_reference_preview end

"""Convert consumed propellant mass into a percent of configured initial mass."""
function rpo_lqmpc_tracking_fuel_used_pct end

"""Track a retimed RPO path with LQ-MPC and report tracking/fuel metrics."""
function rpo_track_retimed_path_lqmpc end

"""Run all configured RPO planners across all comparison cases."""
function rpo_run_planner_comparison_batch(cases, geometry, cfg::RPOPlannerComparisonConfig=RPOPlannerComparisonConfig())
    cfg = _rpo_comparison_config_with_fixed_safe_distance(cfg)
    planners = copy(cfg.planners)
    HYPRSupport.require_hypr()
    # Its Plotly outputs need PlotlyJS: refuse now, not after the batch has run.
    RPOComparisonReports._require_plotly_renderer(cfg.write_plotly_outputs, cfg.write_failed_path_outputs)
    comparison_cases = collect(cases)
    planner_types = [normalize_rpo_comparison_planner_type(p) for p in planners]
    if cfg.optimizer.match_hypr_runtime && (:hypr in planner_types)
        order = vcat([findfirst(==(:hypr), planner_types)], [i for i in eachindex(planner_types) if planner_types[i] != :hypr])
        planners, planner_types = planners[order], planner_types[order]
    end

    results_by_planner = Dict{Symbol, Vector{NamedTuple}}()
    plans_by_planner = Dict{Symbol, Vector{NamedTuple}}()
    master_rng = MersenneTwister(cfg.rng_seed)
    total_runs = length(planner_types) * length(comparison_cases)
    completed_runs = 0
    cfg.show_progress && _rpo_comparison_progress_line!(completed_runs, total_runs)

    for (planner, selection) in zip(planner_types, planners)
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
                selection,
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

include("reporting_compat.jl")

include("executed_tracking.jl")
