"""Grouped PSO settings controlling particle count, iterations, and waypoint search dimensions."""
Base.@kwdef struct RPOPSOSwarmSettings
    n_waypoints::Int = 5
    n_particles::Int = 200
    n_iters::Int = 55
    spread_scale::Float64 = 0.2
    search_margin_m::Float64 = 10.0
    sample_ds_m::Float64 = 0.05
    curve_type::Symbol = :bezier
end

"""Grouped path-cost weights and normalization parameters for RPO HYPR planning."""
Base.@kwdef struct RPOPSOObjectiveSettings
    w_len::Float64 = 1.0
    w_obs::Float64 = 1.0e6
    w_fuel::Float64 = 1.0
    obstacle_sigmoid_k::Float64 = 1.0e6
    obstacle_sigmoid_tol_m::Float64 = 0.0
    w_inertia::Float64 = 0.7
    c1::Float64 = 1.4
    c2::Float64 = 1.4
    cost_ref_distance_m::Float64 = 20.0
    mass_kg::Float64 = 12.0
    tf_s::Float64 = 120.0
    isp_s::Float64 = 60.0
    g0_mps2::Float64 = 9.80665
end

"""Grouped settings for adaptive PSO search-margin and waypoint-count selection."""
Base.@kwdef struct RPOPSOAdaptiveSettings
    enabled::Bool = true
    allow_downscale::Bool = false
    complexity_weight::Float64 = 0.35
    distance_weight::Float64 = 0.65
    waypoint_gain::Float64 = 3.0
    effort_min_fraction::Float64 = 0.75
    effort_max_fraction::Float64 = 1.5
    n_waypoints_min::Int = 3
    n_waypoints_max::Int = 8
    n_particles_min::Int = 140
    n_particles_max::Int = 320
    n_iters_min::Int = 14
    n_iters_max::Int = 55
    w_len_min::Float64 = 0.5
    w_len_max::Float64 = 1.5
    w_obs_min::Float64 = 1.0e6
    w_obs_max::Float64 = 1.0e6
    w_inertia_min::Float64 = 0.4
    w_inertia_max::Float64 = 0.75
    c1_min::Float64 = 1.2
    c1_max::Float64 = 1.8
    c2_min::Float64 = 1.2
    c2_max::Float64 = 2.2
    spread_scale_min::Float64 = 0.05
    spread_scale_max::Float64 = 0.5
end

"""Adaptive collision-sampling controls tied to clearance and curvature."""
Base.@kwdef struct RPOAdaptiveSamplingSettings
    enabled::Bool = true
    max_ds_m::Float64 = 0.50
    far_clearance_m::Float64 = 1.0
    power::Float64 = 1.0
    safe_distance_fraction::Float64 = 0.5
    obstacle_guard_fraction::Float64 = 0.5
end

"""Settings for culling and reseeding weak PSO particles."""
Base.@kwdef struct RPOPSOCullSettings
    enabled::Bool = true
    fraction_max::Float64 = 0.35
    start_iter::Int = 8
    noise_scale::Float64 = 0.25
    arc_velocity_scale::Float64 = 0.12
end

"""Iteration schedule for PSO inertia and acceleration coefficients."""
Base.@kwdef struct RPOPSOScheduleSettings
    enabled::Bool = true
    w_end_fraction::Float64 = 0.65
    c1_end_fraction::Float64 = 0.75
    c2_end_fraction::Float64 = 1.25
    transition_fraction::Float64 = 0.45
    w_min::Float64 = 0.25
    c_min::Float64 = 0.5
    c_max::Float64 = 2.5
end

"""Controls for detecting stagnant PSO progress."""
Base.@kwdef struct RPOPSOStagnationSettings
    stagnation_learning_enable::Bool = true
    stagnation_learning_threshold::Int = 8
    stagnation_learning_elite_fraction::Float64 = 0.10
    stagnation_learning_max_blocks::Int = 2
end

"""Controls for early termination once feasible PSO progress has stabilized."""
Base.@kwdef struct RPOPSOEarlyStoppingSettings
    enabled::Bool = false
    patience::Int = 10
    min_iters::Int = 12
    min_abs_improvement::Float64 = 1.0e-8
    min_rel_improvement::Float64 = 1.0e-4
    require_feasible::Bool = true
end

"""Configuration for quick geometry probes used by adaptive PSO sizing."""
Base.@kwdef struct RPOPSOProbeSettings
    enabled::Bool = true
    max_depth::Int = 2
    candidates::Int = 24
    offset_scale::Float64 = 1.0
    sample_ds_m::Float64 = 0.25
    seed::Int = 1234
end

"""Configuration for widening the search after obstacle-related stagnation."""
Base.@kwdef struct RPOPSOReexploreSettings
    enabled::Bool = true
    trigger_iter::Int = 10
    search_margin_scale::Float64 = 2.0
    waypoint_scale::Float64 = 1.5
    waypoint_increment::Int = 2
    max_waypoints::Int = max(8 + 4, Int(ceil(1.5 * 8)))
end

"""RRT-Connect warm-start settings used before PSO refinement."""
Base.@kwdef struct RPOPSORRTConnectWarmstartSettings
    enabled::Bool = false
    n_iters::Int = 250
    step_size_m::Float64 = 0.75
    goal_sample_rate::Float64 = 0.05
    collision_sample_ds_m::Float64 = 0.10
    connect_max_steps::Int = 10_000
    shortcut_iters::Int = 40
    runtime_limit_s::Float64 = Inf
    box_margin_m::Float64 = 0.75
end

"""Post-PSO shortcutting, handle tightening, and Bezier refit settings."""
Base.@kwdef struct RPOPSORefinementSettings
    enabled::Bool = true
    start_iter::Int = 60
    period::Int = 10
    sample_ds_m::Float64 = 0.025
    merge_distance_m::Float64 = 1.0
    waypoint_passes::Int = 24
    rounds::Int = 16
    min_abs_cost_improvement::Float64 = 1.0e-8
    min_rel_cost_improvement::Float64 = 1.0e-5
    insert_straight_waypoints::Bool = true
    straight_max_segment_length_m::Float64 = 2.0
    straight_max_inserted::Int = 12
    straight_clearance_margin_m::Float64 = 0.0
end

"""Velocity and acceleration limits used when retiming RPO paths."""
Base.@kwdef struct RPOPSORetimingSettings
    dt_s::Float64 = 1.0
    reaction_time_s::Float64 = 0.25
    a_max_mps2::Float64 = 0.02
    speed_scale::Float64 = 0.5
    min_speed_mps::Float64 = 0.0
    max_speed_mps::Float64 = Inf
    max_steps::Int = 100_000
end

"""Hierarchical RPO PSO configuration assembled from grouped setting structs."""
Base.@kwdef struct RPOPSOConfigurator
    swarm::RPOPSOSwarmSettings = RPOPSOSwarmSettings()
    objective::RPOPSOObjectiveSettings = RPOPSOObjectiveSettings()
    adaptive::RPOPSOAdaptiveSettings = RPOPSOAdaptiveSettings()
    adaptive_sampling::RPOAdaptiveSamplingSettings = RPOAdaptiveSamplingSettings()
    cull::RPOPSOCullSettings = RPOPSOCullSettings()
    schedule::RPOPSOScheduleSettings = RPOPSOScheduleSettings()
    stagnation::RPOPSOStagnationSettings = RPOPSOStagnationSettings()
    early_stopping::RPOPSOEarlyStoppingSettings = RPOPSOEarlyStoppingSettings()
    probe::RPOPSOProbeSettings = RPOPSOProbeSettings()
    reexplore::RPOPSOReexploreSettings = RPOPSOReexploreSettings()
    rrt_warmstart::RPOPSORRTConnectWarmstartSettings = RPOPSORRTConnectWarmstartSettings()
    refinement::RPOPSORefinementSettings = RPOPSORefinementSettings()
    retiming::RPOPSORetimingSettings = RPOPSORetimingSettings()
    safe_distance_m::Float64 = 0.0
    goal_collision_margin_m::Float64 = 0.0
end

"""Flattened RPO HYPR/PSO configuration consumed by the planner hot path."""
Base.@kwdef struct RPOPSOConfig
    n_waypoints::Int = 5
    n_particles::Int = 200
    n_iters::Int = 55
    iteration_runtime_limit_s::Float64 = Inf
    w_len::Float64 = 1.0
    w_obs::Float64 = 1.0e6
    w_fuel::Float64 = 1.0
    obstacle_sigmoid_k::Float64 = 1.0e6
    obstacle_sigmoid_tol_m::Float64 = 0.0
    w_inertia::Float64 = 0.7
    c1::Float64 = 1.4
    c2::Float64 = 1.4
    spread_scale::Float64 = 0.2
    search_margin_m::Float64 = 10.0
    sample_ds_m::Float64 = 0.05
    curve_type::Symbol = :bezier
    cost_ref_distance_m::Float64 = 20.0
    mass_kg::Float64 = 12.0
    tf_s::Float64 = 120.0
    isp_s::Float64 = 60.0
    g0_mps2::Float64 = 9.80665
    retime_dt_s::Float64 = 1.0
    retime_reaction_time_s::Float64 = 0.25
    retime_a_max_mps2::Float64 = 0.02
    retime_speed_scale::Float64 = 0.5
    retime_min_speed_mps::Float64 = 0.0
    retime_max_speed_mps::Float64 = Inf
    retime_max_steps::Int = 100_000
    safe_distance_m::Float64 = 0.0
    goal_collision_margin_m::Float64 = 0.0
    adaptive_enable::Bool = true
    adaptive_allow_downscale::Bool = false
    adaptive_complexity_weight::Float64 = 0.35
    adaptive_distance_weight::Float64 = 0.65
    adaptive_waypoint_gain::Float64 = 3.0
    adaptive_effort_min_fraction::Float64 = 0.75
    adaptive_effort_max_fraction::Float64 = 1.5
    adaptive_n_waypoints_min::Int = 3
    adaptive_n_waypoints_max::Int = 8
    adaptive_n_particles_min::Int = 140
    adaptive_n_particles_max::Int = 320
    adaptive_n_iters_min::Int = 14
    adaptive_n_iters_max::Int = 55
    adaptive_w_len_min::Float64 = 0.5
    adaptive_w_len_max::Float64 = 1.5
    adaptive_w_obs_min::Float64 = 1.0e6
    adaptive_w_obs_max::Float64 = 1.0e6
    adaptive_w_inertia_min::Float64 = 0.4
    adaptive_w_inertia_max::Float64 = 0.75
    adaptive_c1_min::Float64 = 1.2
    adaptive_c1_max::Float64 = 1.8
    adaptive_c2_min::Float64 = 1.2
    adaptive_c2_max::Float64 = 2.2
    adaptive_spread_scale_min::Float64 = 0.05
    adaptive_spread_scale_max::Float64 = 0.5
    adaptive_sampling_enable::Bool = true
    adaptive_sampling_max_ds_m::Float64 = 0.50
    adaptive_sampling_far_clearance_m::Float64 = 1.0
    adaptive_sampling_power::Float64 = 1.0
    adaptive_sampling_safe_distance_fraction::Float64 = 0.5
    adaptive_sampling_obstacle_guard_fraction::Float64 = 0.5
    cull_enable::Bool = true
    cull_fraction_max::Float64 = 0.35
    cull_start_iter::Int = 8
    cull_noise_scale::Float64 = 0.25
    cull_arc_velocity_scale::Float64 = 0.12
    schedule_enable::Bool = true
    schedule_w_end_fraction::Float64 = 0.65
    schedule_c1_end_fraction::Float64 = 0.75
    schedule_c2_end_fraction::Float64 = 1.25
    schedule_transition_fraction::Float64 = 0.45
    schedule_w_min::Float64 = 0.25
    schedule_c_min::Float64 = 0.5
    schedule_c_max::Float64 = 2.5
    stagnation_learning_enable::Bool = true
    stagnation_learning_threshold::Int = 8
    stagnation_learning_elite_fraction::Float64 = 0.10
    stagnation_learning_max_blocks::Int = 2
    early_stopping_enable::Bool = false
    early_stopping_patience::Int = 10
    early_stopping_min_iters::Int = 12
    early_stopping_min_abs_improvement::Float64 = 1.0e-8
    early_stopping_min_rel_improvement::Float64 = 1.0e-4
    early_stopping_require_feasible::Bool = true
    probe_enable::Bool = true
    probe_max_depth::Int = 2
    probe_candidates::Int = 24
    probe_offset_scale::Float64 = 1.0
    probe_sample_ds_m::Float64 = 0.25
    probe_seed::Int = 1234
    reexplore_enable::Bool = true
    reexplore_trigger_iter::Int = 10
    reexplore_search_margin_scale::Float64 = 2.0
    reexplore_waypoint_scale::Float64 = 1.5
    reexplore_waypoint_increment::Int = 2
    reexplore_max_waypoints::Int = max(8 + 4, Int(ceil(1.5 * 8)))
    rrt_warmstart_enable::Bool = false
    rrt_warmstart_iters::Int = 250
    rrt_warmstart_step_size_m::Float64 = 0.75
    rrt_warmstart_goal_sample_rate::Float64 = 0.05
    rrt_warmstart_collision_sample_ds_m::Float64 = 0.10
    rrt_warmstart_connect_max_steps::Int = 10_000
    rrt_warmstart_shortcut_iters::Int = 40
    rrt_warmstart_runtime_limit_s::Float64 = Inf
    rrt_warmstart_box_margin_m::Float64 = 0.75
    refinement_enable::Bool = true
    refinement_start_iter::Int = 60
    refinement_period::Int = 10
    refinement_sample_ds_m::Float64 = 0.025
    refinement_merge_distance_m::Float64 = 1.0
    refinement_waypoint_passes::Int = 24
    refinement_rounds::Int = 16
    refinement_min_abs_cost_improvement::Float64 = 1.0e-8
    refinement_min_rel_cost_improvement::Float64 = 1.0e-5
    refinement_insert_straight_waypoints::Bool = true
    refinement_straight_max_segment_length_m::Float64 = 2.0
    refinement_straight_max_inserted::Int = 12
    refinement_straight_clearance_margin_m::Float64 = 0.0
end

"""Convert a flattened RPO PSO config into a named tuple for keyword reconstruction."""
function _rpo_pso_config_tuple(cfg::RPOPSOConfig)
    names = fieldnames(RPOPSOConfig)
    return NamedTuple{names}(map(name -> getfield(cfg, name), names))
end

"""Fill sampling-distance defaults from the configured keepout distance when needed."""
function _rpo_pso_sync_sample_ds_with_safe_distance(values)
    safe_distance_m = Float64(get(values, :safe_distance_m, 0.0))
    safe_distance_m > 0.0 || return values
    return merge(values, (sample_ds_m=safe_distance_m,))
end

"""Return the collision-sampling spacing used for primary RPO HYPR evaluation."""
function rpo_hypr_sampling_density_m(cfg::RPOPSOConfig, safe_distance_m::Real=cfg.safe_distance_m)
    safe_distance = Float64(safe_distance_m)
    safe_distance > 0.0 && return safe_distance
    cfg.safe_distance_m > 0.0 && return cfg.safe_distance_m
    return cfg.sample_ds_m
end

"""Return the collision-sampling spacing used during RPO post-refinement."""
function rpo_hypr_refinement_sampling_density_m(cfg::RPOPSOConfig, safe_distance_m::Real=cfg.safe_distance_m)
    safe_distance = Float64(safe_distance_m)
    safe_distance > 0.0 && return safe_distance
    cfg.safe_distance_m > 0.0 && return cfg.safe_distance_m
    return cfg.refinement_sample_ds_m
end

const RPO_PSO_CONFIG_ALIASES = Dict{Symbol, Symbol}(
    :pso_waypoints => :n_waypoints,
    :pso_particles => :n_particles,
    :pso_iters => :n_iters,
    :pso_iteration_runtime_limit => :iteration_runtime_limit_s,
    :pso_iteration_runtime_limit_s => :iteration_runtime_limit_s,
    :pso_w_len => :w_len,
    :pso_w_obs => :w_obs,
    :pso_w_fuel => :w_fuel,
    :pso_obstacle_sigmoid_k => :obstacle_sigmoid_k,
    :pso_obstacle_sigmoid_tol => :obstacle_sigmoid_tol_m,
    :pso_obstacle_sigmoid_tol_m => :obstacle_sigmoid_tol_m,
    :pso_w_inertia => :w_inertia,
    :pso_c1 => :c1,
    :pso_c2 => :c2,
    :pso_cull_fraction_max => :cull_fraction_max,
    :pso_cull_start_iter => :cull_start_iter,
    :pso_cull_noise_scale => :cull_noise_scale,
    :pso_cull_arc_velocity_scale => :cull_arc_velocity_scale,
    :pso_schedule_enable => :schedule_enable,
    :pso_schedule_w_end_fraction => :schedule_w_end_fraction,
    :pso_schedule_c1_end_fraction => :schedule_c1_end_fraction,
    :pso_schedule_c2_end_fraction => :schedule_c2_end_fraction,
    :pso_schedule_transition_fraction => :schedule_transition_fraction,
    :pso_schedule_w_min => :schedule_w_min,
    :pso_schedule_c_min => :schedule_c_min,
    :pso_schedule_c_max => :schedule_c_max,
    :pso_stagnation_learning_enable => :stagnation_learning_enable,
    :pso_stagnation_learning_threshold => :stagnation_learning_threshold,
    :pso_stagnation_learning_elite_fraction => :stagnation_learning_elite_fraction,
    :pso_stagnation_learning_max_blocks => :stagnation_learning_max_blocks,
    :pso_early_stopping_enable => :early_stopping_enable,
    :pso_early_stopping_patience => :early_stopping_patience,
    :pso_early_stopping_min_iters => :early_stopping_min_iters,
    :pso_early_stopping_min_abs_improvement => :early_stopping_min_abs_improvement,
    :pso_early_stopping_min_rel_improvement => :early_stopping_min_rel_improvement,
    :pso_early_stopping_require_feasible => :early_stopping_require_feasible,
    :pso_safe_distance => :safe_distance_m,
    :pso_goal_collision_margin => :goal_collision_margin_m,
    :pso_search_margin => :search_margin_m,
    :pso_spread_scale => :spread_scale,
    :pso_spread_scale_min => :adaptive_spread_scale_min,
    :pso_spread_scale_max => :adaptive_spread_scale_max,
    :pso_sample_ds => :sample_ds_m,
    :pso_curve_type => :curve_type,
    :pso_adaptive => :adaptive_enable,
    :pso_adaptive_allow_downscale => :adaptive_allow_downscale,
    :pso_adaptive_sampling => :adaptive_sampling_enable,
    :pso_adaptive_sampling_enable => :adaptive_sampling_enable,
    :pso_adaptive_sampling_max_ds => :adaptive_sampling_max_ds_m,
    :pso_adaptive_sampling_max_ds_m => :adaptive_sampling_max_ds_m,
    :pso_adaptive_sampling_far_clearance => :adaptive_sampling_far_clearance_m,
    :pso_adaptive_sampling_far_clearance_m => :adaptive_sampling_far_clearance_m,
    :pso_adaptive_sampling_power => :adaptive_sampling_power,
    :pso_adaptive_sampling_safe_distance_fraction => :adaptive_sampling_safe_distance_fraction,
    :pso_adaptive_sampling_obstacle_guard_fraction => :adaptive_sampling_obstacle_guard_fraction,
    :pso_probe_enable => :probe_enable,
    :pso_probe_max_depth => :probe_max_depth,
    :pso_probe_candidates => :probe_candidates,
    :pso_probe_offset_scale => :probe_offset_scale,
    :pso_probe_sample_ds => :probe_sample_ds_m,
    :pso_probe_seed => :probe_seed,
    :pso_cost_ref_distance => :cost_ref_distance_m,
    :pso_complexity_weight => :adaptive_complexity_weight,
    :pso_w_inertia_min => :adaptive_w_inertia_min,
    :pso_w_inertia_max => :adaptive_w_inertia_max,
    :pso_c1_min => :adaptive_c1_min,
    :pso_c1_max => :adaptive_c1_max,
    :pso_c2_min => :adaptive_c2_min,
    :pso_c2_max => :adaptive_c2_max,
    :pso_w_obs_min => :adaptive_w_obs_min,
    :pso_w_obs_max => :adaptive_w_obs_max,
    :pso_w_len_min => :adaptive_w_len_min,
    :pso_w_len_max => :adaptive_w_len_max,
    :pso_particles_min => :adaptive_n_particles_min,
    :pso_particles_max => :adaptive_n_particles_max,
    :pso_iters_min => :adaptive_n_iters_min,
    :pso_iters_max => :adaptive_n_iters_max,
    :pso_waypoints_min => :adaptive_n_waypoints_min,
    :pso_waypoints_max => :adaptive_n_waypoints_max,
    :pso_reexplore_enable => :reexplore_enable,
    :pso_reexplore_trigger_iter => :reexplore_trigger_iter,
    :pso_reexplore_search_margin_scale => :reexplore_search_margin_scale,
    :pso_reexplore_waypoint_scale => :reexplore_waypoint_scale,
    :pso_reexplore_waypoint_increment => :reexplore_waypoint_increment,
    :pso_reexplore_max_waypoints => :reexplore_max_waypoints,
    :pso_rrt_warmstart_enable => :rrt_warmstart_enable,
    :pso_rrt_warmstart_iters => :rrt_warmstart_iters,
    :pso_rrt_warmstart_step_size => :rrt_warmstart_step_size_m,
    :pso_rrt_warmstart_goal_sample_rate => :rrt_warmstart_goal_sample_rate,
    :pso_rrt_warmstart_collision_sample_ds => :rrt_warmstart_collision_sample_ds_m,
    :pso_rrt_warmstart_connect_max_steps => :rrt_warmstart_connect_max_steps,
    :pso_rrt_warmstart_shortcut_iters => :rrt_warmstart_shortcut_iters,
    :pso_rrt_warmstart_runtime_limit => :rrt_warmstart_runtime_limit_s,
    :pso_rrt_warmstart_box_margin => :rrt_warmstart_box_margin_m,
    :pso_rrt_warmstart_box_margin_m => :rrt_warmstart_box_margin_m,
)

"""Apply supported RPO PSO keyword aliases before constructing the flattened config."""
function _rpo_pso_normalize_kwargs(kwargs)
    pairs = Pair{Symbol, Any}[]
    for (key, value) in kwargs
        push!(pairs, get(RPO_PSO_CONFIG_ALIASES, key, key) => value)
    end
    return (; pairs...)
end

"""Flattened RPO HYPR/PSO configuration consumed by the planner hot path."""
function RPOPSOConfig(configurator::RPOPSOConfigurator; kwargs...)
    swarm = configurator.swarm
    objective = configurator.objective
    adaptive = configurator.adaptive
    adaptive_sampling = configurator.adaptive_sampling
    cull = configurator.cull
    schedule = configurator.schedule
    stagnation = configurator.stagnation
    early_stopping = configurator.early_stopping
    probe = configurator.probe
    reexplore = configurator.reexplore
    rrt_warmstart = configurator.rrt_warmstart
    refinement = configurator.refinement
    retiming = configurator.retiming
    cfg = RPOPSOConfig(;
        n_waypoints=swarm.n_waypoints,
        n_particles=swarm.n_particles,
        n_iters=swarm.n_iters,
        w_len=objective.w_len,
        w_obs=objective.w_obs,
        w_fuel=objective.w_fuel,
        obstacle_sigmoid_k=objective.obstacle_sigmoid_k,
        obstacle_sigmoid_tol_m=objective.obstacle_sigmoid_tol_m,
        w_inertia=objective.w_inertia,
        c1=objective.c1,
        c2=objective.c2,
        spread_scale=swarm.spread_scale,
        search_margin_m=swarm.search_margin_m,
        sample_ds_m=swarm.sample_ds_m,
        curve_type=swarm.curve_type,
        cost_ref_distance_m=objective.cost_ref_distance_m,
        mass_kg=objective.mass_kg,
        tf_s=objective.tf_s,
        isp_s=objective.isp_s,
        g0_mps2=objective.g0_mps2,
        retime_dt_s=retiming.dt_s,
        retime_reaction_time_s=retiming.reaction_time_s,
        retime_a_max_mps2=retiming.a_max_mps2,
        retime_speed_scale=retiming.speed_scale,
        retime_min_speed_mps=retiming.min_speed_mps,
        retime_max_speed_mps=retiming.max_speed_mps,
        retime_max_steps=retiming.max_steps,
        safe_distance_m=configurator.safe_distance_m,
        goal_collision_margin_m=configurator.goal_collision_margin_m,
        adaptive_enable=adaptive.enabled,
        adaptive_allow_downscale=adaptive.allow_downscale,
        adaptive_complexity_weight=adaptive.complexity_weight,
        adaptive_distance_weight=adaptive.distance_weight,
        adaptive_waypoint_gain=adaptive.waypoint_gain,
        adaptive_effort_min_fraction=adaptive.effort_min_fraction,
        adaptive_effort_max_fraction=adaptive.effort_max_fraction,
        adaptive_n_waypoints_min=adaptive.n_waypoints_min,
        adaptive_n_waypoints_max=adaptive.n_waypoints_max,
        adaptive_n_particles_min=adaptive.n_particles_min,
        adaptive_n_particles_max=adaptive.n_particles_max,
        adaptive_n_iters_min=adaptive.n_iters_min,
        adaptive_n_iters_max=adaptive.n_iters_max,
        adaptive_w_len_min=adaptive.w_len_min,
        adaptive_w_len_max=adaptive.w_len_max,
        adaptive_w_obs_min=adaptive.w_obs_min,
        adaptive_w_obs_max=adaptive.w_obs_max,
        adaptive_w_inertia_min=adaptive.w_inertia_min,
        adaptive_w_inertia_max=adaptive.w_inertia_max,
        adaptive_c1_min=adaptive.c1_min,
        adaptive_c1_max=adaptive.c1_max,
        adaptive_c2_min=adaptive.c2_min,
        adaptive_c2_max=adaptive.c2_max,
        adaptive_spread_scale_min=adaptive.spread_scale_min,
        adaptive_spread_scale_max=adaptive.spread_scale_max,
        adaptive_sampling_enable=adaptive_sampling.enabled,
        adaptive_sampling_max_ds_m=adaptive_sampling.max_ds_m,
        adaptive_sampling_far_clearance_m=adaptive_sampling.far_clearance_m,
        adaptive_sampling_power=adaptive_sampling.power,
        adaptive_sampling_safe_distance_fraction=adaptive_sampling.safe_distance_fraction,
        adaptive_sampling_obstacle_guard_fraction=adaptive_sampling.obstacle_guard_fraction,
        cull_enable=cull.enabled,
        cull_fraction_max=cull.fraction_max,
        cull_start_iter=cull.start_iter,
        cull_noise_scale=cull.noise_scale,
        cull_arc_velocity_scale=cull.arc_velocity_scale,
        schedule_enable=schedule.enabled,
        schedule_w_end_fraction=schedule.w_end_fraction,
        schedule_c1_end_fraction=schedule.c1_end_fraction,
        schedule_c2_end_fraction=schedule.c2_end_fraction,
        schedule_transition_fraction=schedule.transition_fraction,
        schedule_w_min=schedule.w_min,
        schedule_c_min=schedule.c_min,
        schedule_c_max=schedule.c_max,
        stagnation_learning_enable=stagnation.stagnation_learning_enable,
        stagnation_learning_threshold=stagnation.stagnation_learning_threshold,
        stagnation_learning_elite_fraction=stagnation.stagnation_learning_elite_fraction,
        stagnation_learning_max_blocks=stagnation.stagnation_learning_max_blocks,
        early_stopping_enable=early_stopping.enabled,
        early_stopping_patience=early_stopping.patience,
        early_stopping_min_iters=early_stopping.min_iters,
        early_stopping_min_abs_improvement=early_stopping.min_abs_improvement,
        early_stopping_min_rel_improvement=early_stopping.min_rel_improvement,
        early_stopping_require_feasible=early_stopping.require_feasible,
        probe_enable=probe.enabled,
        probe_max_depth=probe.max_depth,
        probe_candidates=probe.candidates,
        probe_offset_scale=probe.offset_scale,
        probe_sample_ds_m=probe.sample_ds_m,
        probe_seed=probe.seed,
        reexplore_enable=reexplore.enabled,
        reexplore_trigger_iter=reexplore.trigger_iter,
        reexplore_search_margin_scale=reexplore.search_margin_scale,
        reexplore_waypoint_scale=reexplore.waypoint_scale,
        reexplore_waypoint_increment=reexplore.waypoint_increment,
        reexplore_max_waypoints=reexplore.max_waypoints,
        rrt_warmstart_enable=rrt_warmstart.enabled,
        rrt_warmstart_iters=rrt_warmstart.n_iters,
        rrt_warmstart_step_size_m=rrt_warmstart.step_size_m,
        rrt_warmstart_goal_sample_rate=rrt_warmstart.goal_sample_rate,
        rrt_warmstart_collision_sample_ds_m=rrt_warmstart.collision_sample_ds_m,
        rrt_warmstart_connect_max_steps=rrt_warmstart.connect_max_steps,
        rrt_warmstart_shortcut_iters=rrt_warmstart.shortcut_iters,
        rrt_warmstart_runtime_limit_s=rrt_warmstart.runtime_limit_s,
        rrt_warmstart_box_margin_m=rrt_warmstart.box_margin_m,
        refinement_enable=refinement.enabled,
        refinement_start_iter=refinement.start_iter,
        refinement_period=refinement.period,
        refinement_sample_ds_m=refinement.sample_ds_m,
        refinement_merge_distance_m=refinement.merge_distance_m,
        refinement_waypoint_passes=refinement.waypoint_passes,
        refinement_rounds=refinement.rounds,
        refinement_min_abs_cost_improvement=refinement.min_abs_cost_improvement,
        refinement_min_rel_cost_improvement=refinement.min_rel_cost_improvement,
        refinement_insert_straight_waypoints=refinement.insert_straight_waypoints,
        refinement_straight_max_segment_length_m=refinement.straight_max_segment_length_m,
        refinement_straight_max_inserted=refinement.straight_max_inserted,
        refinement_straight_clearance_margin_m=refinement.straight_clearance_margin_m,
    )
    return rpo_pso_config(cfg; kwargs...)
end

"""Construct an RPO PSO config by overriding an existing config or configurator."""
function rpo_pso_config(base::RPOPSOConfig=RPOPSOConfig(); kwargs...)
    values = merge(_rpo_pso_config_tuple(base), _rpo_pso_normalize_kwargs(kwargs))
    values = _rpo_pso_sync_sample_ds_with_safe_distance(values)
    return validate_rpo_pso_config(RPOPSOConfig(; values...))
end

"""Construct an RPO PSO config by overriding an existing config or configurator."""
rpo_pso_config(configurator::RPOPSOConfigurator; kwargs...) = RPOPSOConfig(configurator; kwargs...)

"""Validate RPO PSO configuration ranges before planner use."""
function validate_rpo_pso_config(cfg::RPOPSOConfig)
    cfg.n_waypoints >= 0 || throw(ArgumentError("n_waypoints must be nonnegative."))
    cfg.n_particles > 0 || throw(ArgumentError("n_particles must be positive."))
    cfg.n_iters >= 0 || throw(ArgumentError("n_iters must be nonnegative."))
    cfg.iteration_runtime_limit_s >= 0.0 ||
        throw(ArgumentError("iteration_runtime_limit_s must be nonnegative."))
    cfg.sample_ds_m > 0.0 || throw(ArgumentError("sample_ds_m must be positive."))
    cfg.curve_type in (:bezier, :polyline) ||
        throw(ArgumentError("curve_type must be :bezier or :polyline."))
    cfg.obstacle_sigmoid_k > 0.0 || throw(ArgumentError("obstacle_sigmoid_k must be positive."))
    cfg.obstacle_sigmoid_tol_m >= 0.0 || throw(ArgumentError("obstacle_sigmoid_tol_m must be nonnegative."))
    cfg.cost_ref_distance_m >= 0.0 || throw(ArgumentError("cost_ref_distance_m must be nonnegative."))
    cfg.mass_kg > 0.0 || throw(ArgumentError("mass_kg must be positive."))
    cfg.tf_s > 0.0 || throw(ArgumentError("tf_s must be positive."))
    cfg.isp_s > 0.0 || throw(ArgumentError("isp_s must be positive."))
    cfg.g0_mps2 > 0.0 || throw(ArgumentError("g0_mps2 must be positive."))
    cfg.retime_dt_s > 0.0 || throw(ArgumentError("retime_dt_s must be positive."))
    cfg.retime_a_max_mps2 > 0.0 || throw(ArgumentError("retime_a_max_mps2 must be positive."))
    cfg.retime_speed_scale > 0.0 || throw(ArgumentError("retime_speed_scale must be positive."))
    cfg.retime_min_speed_mps >= 0.0 || throw(ArgumentError("retime_min_speed_mps must be nonnegative."))
    cfg.retime_max_speed_mps >= cfg.retime_min_speed_mps ||
        throw(ArgumentError("retime_max_speed_mps must be at least retime_min_speed_mps."))
    cfg.retime_max_steps > 0 || throw(ArgumentError("retime_max_steps must be positive."))
    cfg.safe_distance_m >= 0.0 || throw(ArgumentError("safe_distance_m must be nonnegative."))
    cfg.goal_collision_margin_m >= 0.0 || throw(ArgumentError("goal_collision_margin_m must be nonnegative."))
    0.0 <= cfg.adaptive_complexity_weight || throw(ArgumentError("adaptive_complexity_weight must be nonnegative."))
    0.0 <= cfg.adaptive_distance_weight || throw(ArgumentError("adaptive_distance_weight must be nonnegative."))
    cfg.adaptive_effort_min_fraction > 0.0 || throw(ArgumentError("adaptive_effort_min_fraction must be positive."))
    cfg.adaptive_effort_max_fraction >= cfg.adaptive_effort_min_fraction ||
        throw(ArgumentError("adaptive_effort_max_fraction must be at least adaptive_effort_min_fraction."))
    cfg.adaptive_n_waypoints_min >= 0 || throw(ArgumentError("adaptive_n_waypoints_min must be nonnegative."))
    cfg.adaptive_n_waypoints_max >= 0 || throw(ArgumentError("adaptive_n_waypoints_max must be nonnegative."))
    cfg.adaptive_n_particles_min > 0 || throw(ArgumentError("adaptive_n_particles_min must be positive."))
    cfg.adaptive_n_particles_max >= 0 || throw(ArgumentError("adaptive_n_particles_max must be nonnegative."))
    cfg.adaptive_n_iters_min >= 0 || throw(ArgumentError("adaptive_n_iters_min must be nonnegative."))
    cfg.adaptive_n_iters_max >= 0 || throw(ArgumentError("adaptive_n_iters_max must be nonnegative."))
    cfg.adaptive_w_len_min <= cfg.adaptive_w_len_max ||
        throw(ArgumentError("adaptive_w_len_min must be at most adaptive_w_len_max."))
    cfg.adaptive_w_obs_min <= cfg.adaptive_w_obs_max ||
        throw(ArgumentError("adaptive_w_obs_min must be at most adaptive_w_obs_max."))
    cfg.adaptive_w_inertia_min <= cfg.adaptive_w_inertia_max ||
        throw(ArgumentError("adaptive_w_inertia_min must be at most adaptive_w_inertia_max."))
    cfg.adaptive_c1_min <= cfg.adaptive_c1_max ||
        throw(ArgumentError("adaptive_c1_min must be at most adaptive_c1_max."))
    cfg.adaptive_c2_min <= cfg.adaptive_c2_max ||
        throw(ArgumentError("adaptive_c2_min must be at most adaptive_c2_max."))
    cfg.adaptive_spread_scale_min <= cfg.adaptive_spread_scale_max ||
        throw(ArgumentError("adaptive_spread_scale_min must be at most adaptive_spread_scale_max."))
    cfg.adaptive_sampling_max_ds_m > 0.0 ||
        throw(ArgumentError("adaptive_sampling_max_ds_m must be positive."))
    cfg.adaptive_sampling_far_clearance_m > 0.0 ||
        throw(ArgumentError("adaptive_sampling_far_clearance_m must be positive."))
    cfg.adaptive_sampling_power > 0.0 ||
        throw(ArgumentError("adaptive_sampling_power must be positive."))
    cfg.adaptive_sampling_safe_distance_fraction > 0.0 ||
        throw(ArgumentError("adaptive_sampling_safe_distance_fraction must be positive."))
    0.0 < cfg.adaptive_sampling_obstacle_guard_fraction <= 1.0 ||
        throw(ArgumentError("adaptive_sampling_obstacle_guard_fraction must be in (0, 1]."))
    0.0 <= cfg.cull_fraction_max <= 1.0 ||
        throw(ArgumentError("cull_fraction_max must be between 0 and 1."))
    cfg.cull_start_iter >= 0 || throw(ArgumentError("cull_start_iter must be nonnegative."))
    cfg.cull_noise_scale >= 0.0 || throw(ArgumentError("cull_noise_scale must be nonnegative."))
    cfg.cull_arc_velocity_scale >= 0.0 || throw(ArgumentError("cull_arc_velocity_scale must be nonnegative."))
    cfg.schedule_transition_fraction > 0.0 || throw(ArgumentError("schedule_transition_fraction must be positive."))
    cfg.schedule_w_min >= 0.0 || throw(ArgumentError("schedule_w_min must be nonnegative."))
    cfg.schedule_c_min >= 0.0 || throw(ArgumentError("schedule_c_min must be nonnegative."))
    cfg.schedule_c_min <= cfg.schedule_c_max ||
        throw(ArgumentError("schedule_c_min must be at most schedule_c_max."))
    cfg.stagnation_learning_threshold >= 0 ||
        throw(ArgumentError("stagnation_learning_threshold must be nonnegative."))
    0.0 <= cfg.stagnation_learning_elite_fraction <= 1.0 ||
        throw(ArgumentError("stagnation_learning_elite_fraction must be between 0 and 1."))
    cfg.stagnation_learning_max_blocks >= 0 ||
        throw(ArgumentError("stagnation_learning_max_blocks must be nonnegative."))
    cfg.early_stopping_patience >= 0 || throw(ArgumentError("early_stopping_patience must be nonnegative."))
    cfg.early_stopping_min_iters >= 0 || throw(ArgumentError("early_stopping_min_iters must be nonnegative."))
    cfg.early_stopping_min_abs_improvement >= 0.0 ||
        throw(ArgumentError("early_stopping_min_abs_improvement must be nonnegative."))
    cfg.early_stopping_min_rel_improvement >= 0.0 ||
        throw(ArgumentError("early_stopping_min_rel_improvement must be nonnegative."))
    cfg.probe_max_depth >= 0 || throw(ArgumentError("probe_max_depth must be nonnegative."))
    cfg.probe_candidates >= 0 || throw(ArgumentError("probe_candidates must be nonnegative."))
    cfg.probe_offset_scale >= 0.0 || throw(ArgumentError("probe_offset_scale must be nonnegative."))
    cfg.probe_sample_ds_m > 0.0 || throw(ArgumentError("probe_sample_ds_m must be positive."))
    cfg.reexplore_trigger_iter >= 0 || throw(ArgumentError("reexplore_trigger_iter must be nonnegative."))
    cfg.reexplore_search_margin_scale > 0.0 ||
        throw(ArgumentError("reexplore_search_margin_scale must be positive."))
    cfg.reexplore_waypoint_scale > 0.0 ||
        throw(ArgumentError("reexplore_waypoint_scale must be positive."))
    cfg.reexplore_waypoint_increment >= 0 ||
        throw(ArgumentError("reexplore_waypoint_increment must be nonnegative."))
    cfg.reexplore_max_waypoints >= 0 || throw(ArgumentError("reexplore_max_waypoints must be nonnegative."))
    cfg.rrt_warmstart_iters >= 0 || throw(ArgumentError("rrt_warmstart_iters must be nonnegative."))
    cfg.rrt_warmstart_step_size_m > 0.0 ||
        throw(ArgumentError("rrt_warmstart_step_size_m must be positive."))
    0.0 <= cfg.rrt_warmstart_goal_sample_rate <= 1.0 ||
        throw(ArgumentError("rrt_warmstart_goal_sample_rate must be between 0 and 1."))
    cfg.rrt_warmstart_collision_sample_ds_m > 0.0 ||
        throw(ArgumentError("rrt_warmstart_collision_sample_ds_m must be positive."))
    cfg.rrt_warmstart_connect_max_steps > 0 ||
        throw(ArgumentError("rrt_warmstart_connect_max_steps must be positive."))
    cfg.rrt_warmstart_shortcut_iters >= 0 ||
        throw(ArgumentError("rrt_warmstart_shortcut_iters must be nonnegative."))
    cfg.rrt_warmstart_runtime_limit_s >= 0.0 ||
        throw(ArgumentError("rrt_warmstart_runtime_limit_s must be nonnegative."))
    cfg.rrt_warmstart_box_margin_m >= 0.0 ||
        throw(ArgumentError("rrt_warmstart_box_margin_m must be nonnegative."))
    cfg.refinement_start_iter >= 0 || throw(ArgumentError("refinement_start_iter must be nonnegative."))
    cfg.refinement_period > 0 || throw(ArgumentError("refinement_period must be positive."))
    cfg.refinement_sample_ds_m > 0.0 || throw(ArgumentError("refinement_sample_ds_m must be positive."))
    cfg.refinement_merge_distance_m >= 0.0 ||
        throw(ArgumentError("refinement_merge_distance_m must be nonnegative."))
    cfg.refinement_waypoint_passes >= 0 ||
        throw(ArgumentError("refinement_waypoint_passes must be nonnegative."))
    cfg.refinement_rounds >= 0 || throw(ArgumentError("refinement_rounds must be nonnegative."))
    cfg.refinement_min_abs_cost_improvement >= 0.0 ||
        throw(ArgumentError("refinement_min_abs_cost_improvement must be nonnegative."))
    cfg.refinement_min_rel_cost_improvement >= 0.0 ||
        throw(ArgumentError("refinement_min_rel_cost_improvement must be nonnegative."))
    cfg.refinement_straight_max_segment_length_m > 0.0 ||
        throw(ArgumentError("refinement_straight_max_segment_length_m must be positive."))
    cfg.refinement_straight_max_inserted >= 0 ||
        throw(ArgumentError("refinement_straight_max_inserted must be nonnegative."))
    cfg.refinement_straight_clearance_margin_m >= 0.0 ||
        throw(ArgumentError("refinement_straight_clearance_margin_m must be nonnegative."))
    return cfg
end
