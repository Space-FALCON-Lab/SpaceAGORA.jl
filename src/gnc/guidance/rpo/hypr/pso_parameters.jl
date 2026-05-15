Base.@kwdef struct RPOPSOSwarmSettings
    n_waypoints::Int = 5
    n_particles::Int = 200
    n_iters::Int = 55
    spread_scale::Float64 = 0.2
    search_margin_m::Float64 = 6.0
    sample_ds_m::Float64 = 0.05
    curve_type::Symbol = :bezier
end

Base.@kwdef struct RPOPSOObjectiveSettings
    w_len::Float64 = 1.0
    w_detour::Float64 = 0.0
    w_dir::Float64 = 0.0
    w_obs::Float64 = 1.0e6
    w_fuel::Float64 = 1.0
    fuel_proxy_mode::Symbol = :retimed
    w_inertia::Float64 = 0.7
    c1::Float64 = 1.4
    c2::Float64 = 1.4
    detour_ratio_ref::Float64 = 2.0
    dist_ref_m::Float64 = 20.0
    cost_ref_distance_m::Float64 = 20.0
    mass_kg::Float64 = 12.0
    tf_s::Float64 = 120.0
    isp_s::Float64 = 60.0
    g0_mps2::Float64 = 9.80665
end

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

Base.@kwdef struct RPOPSOCullSettings
    enabled::Bool = true
    fraction_max::Float64 = 0.35
    start_iter::Int = 8
    noise_scale::Float64 = 0.25
    corner_pull::Float64 = 0.18
    obstacle_pull::Float64 = 0.12
    arc_velocity_scale::Float64 = 0.12
end

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

Base.@kwdef struct RPOPSOStagnationSettings
    stagnation_learning_enable::Bool = true
    stagnation_learning_threshold::Int = 8
    stagnation_learning_elite_fraction::Float64 = 0.10
    stagnation_learning_max_blocks::Int = 2
end

Base.@kwdef struct RPOPSOProbeSettings
    enabled::Bool = true
    max_depth::Int = 2
    candidates::Int = 24
    offset_scale::Float64 = 1.0
    sample_ds_m::Float64 = 0.25
    seed::Int = 1234
end

Base.@kwdef struct RPOPSOReexploreSettings
    enabled::Bool = true
    trigger_iter::Int = 10
    search_margin_scale::Float64 = 1.5
    waypoint_scale::Float64 = 1.5
    waypoint_increment::Int = 2
    max_waypoints::Int = max(8 + 4, Int(ceil(1.5 * 8)))
end

Base.@kwdef struct RPOPSOGIFSettings
    single_run_enable::Bool = true
    top_k::Int = 10
    fps::Int = 6
end

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

Base.@kwdef struct RPOPSORetimingSettings
    dt_s::Float64 = 1.0
    reaction_time_s::Float64 = 0.25
    a_max_mps2::Float64 = 0.02
    speed_scale::Float64 = 0.5
    min_speed_mps::Float64 = 0.0
    max_speed_mps::Float64 = Inf
    max_steps::Int = 100_000
end

Base.@kwdef struct RPOPSOConfigurator
    swarm::RPOPSOSwarmSettings = RPOPSOSwarmSettings()
    objective::RPOPSOObjectiveSettings = RPOPSOObjectiveSettings()
    adaptive::RPOPSOAdaptiveSettings = RPOPSOAdaptiveSettings()
    cull::RPOPSOCullSettings = RPOPSOCullSettings()
    schedule::RPOPSOScheduleSettings = RPOPSOScheduleSettings()
    stagnation::RPOPSOStagnationSettings = RPOPSOStagnationSettings()
    probe::RPOPSOProbeSettings = RPOPSOProbeSettings()
    reexplore::RPOPSOReexploreSettings = RPOPSOReexploreSettings()
    gif::RPOPSOGIFSettings = RPOPSOGIFSettings()
    refinement::RPOPSORefinementSettings = RPOPSORefinementSettings()
    retiming::RPOPSORetimingSettings = RPOPSORetimingSettings()
    safe_distance_m::Float64 = 0.0
    goal_collision_margin_m::Float64 = 0.0
end

Base.@kwdef struct RPOPSOConfig
    n_waypoints::Int = 5
    n_particles::Int = 200
    n_iters::Int = 55
    w_len::Float64 = 1.0
    w_detour::Float64 = 0.0
    w_dir::Float64 = 0.0
    w_obs::Float64 = 1.0e6
    w_fuel::Float64 = 1.0
    fuel_proxy_mode::Symbol = :retimed
    w_inertia::Float64 = 0.7
    c1::Float64 = 1.4
    c2::Float64 = 1.4
    spread_scale::Float64 = 0.2
    search_margin_m::Float64 = 6.0
    sample_ds_m::Float64 = 0.05
    curve_type::Symbol = :bezier
    detour_ratio_ref::Float64 = 2.0
    dist_ref_m::Float64 = 20.0
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
    cull_enable::Bool = true
    cull_fraction_max::Float64 = 0.35
    cull_start_iter::Int = 8
    cull_noise_scale::Float64 = 0.25
    cull_corner_pull::Float64 = 0.18
    cull_obstacle_pull::Float64 = 0.12
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
    probe_enable::Bool = true
    probe_max_depth::Int = 2
    probe_candidates::Int = 24
    probe_offset_scale::Float64 = 1.0
    probe_sample_ds_m::Float64 = 0.25
    probe_seed::Int = 1234
    reexplore_enable::Bool = true
    reexplore_trigger_iter::Int = 10
    reexplore_search_margin_scale::Float64 = 1.5
    reexplore_waypoint_scale::Float64 = 1.5
    reexplore_waypoint_increment::Int = 2
    reexplore_max_waypoints::Int = max(8 + 4, Int(ceil(1.5 * 8)))
    single_run_gif_enable::Bool = true
    gif_top_k::Int = 10
    gif_fps::Int = 6
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

function _rpo_pso_config_tuple(cfg::RPOPSOConfig)
    names = fieldnames(RPOPSOConfig)
    return NamedTuple{names}(map(name -> getfield(cfg, name), names))
end

const RPO_PSO_CONFIG_ALIASES = Dict{Symbol, Symbol}(
    :pso_waypoints => :n_waypoints,
    :pso_particles => :n_particles,
    :pso_iters => :n_iters,
    :pso_w_len => :w_len,
    :pso_w_detour => :w_detour,
    :pso_w_dir => :w_dir,
    :pso_w_obs => :w_obs,
    :pso_w_fuel => :w_fuel,
    :pso_fuel_proxy_mode => :fuel_proxy_mode,
    :pso_w_inertia => :w_inertia,
    :pso_c1 => :c1,
    :pso_c2 => :c2,
    :pso_cull_fraction_max => :cull_fraction_max,
    :pso_cull_start_iter => :cull_start_iter,
    :pso_cull_noise_scale => :cull_noise_scale,
    :pso_cull_corner_pull => :cull_corner_pull,
    :pso_cull_obstacle_pull => :cull_obstacle_pull,
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
    :pso_safe_distance => :safe_distance_m,
    :pso_goal_collision_margin => :goal_collision_margin_m,
    :pso_search_margin => :search_margin_m,
    :pso_spread_scale => :spread_scale,
    :pso_spread_scale_min => :adaptive_spread_scale_min,
    :pso_spread_scale_max => :adaptive_spread_scale_max,
    :pso_sample_ds => :sample_ds_m,
    :pso_curve_type => :curve_type,
    :pso_adaptive => :adaptive_enable,
    :pso_detour_ratio_ref => :detour_ratio_ref,
    :pso_probe_enable => :probe_enable,
    :pso_probe_max_depth => :probe_max_depth,
    :pso_probe_candidates => :probe_candidates,
    :pso_probe_offset_scale => :probe_offset_scale,
    :pso_probe_sample_ds => :probe_sample_ds_m,
    :pso_probe_seed => :probe_seed,
    :pso_dist_ref => :dist_ref_m,
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
    :pso_single_run_gif_enable => :single_run_gif_enable,
    :pso_gif_top_k => :gif_top_k,
    :pso_gif_fps => :gif_fps,
    :pso_shortcut_enable => :refinement_enable,
    :pso_shortcut_start_iter => :refinement_start_iter,
    :pso_shortcut_period => :refinement_period,
    :pso_shortcut_sample_ds => :refinement_sample_ds_m,
    :pso_shortcut_merge_distance => :refinement_merge_distance_m,
    :pso_shortcut_waypoint_passes => :refinement_waypoint_passes,
    :pso_shortcut_rounds => :refinement_rounds,
    :pso_shortcut_min_abs_cost_improvement => :refinement_min_abs_cost_improvement,
    :pso_shortcut_min_rel_cost_improvement => :refinement_min_rel_cost_improvement,
    :pso_shortcut_insert_straight_waypoints => :refinement_insert_straight_waypoints,
    :pso_shortcut_straight_max_segment_length => :refinement_straight_max_segment_length_m,
    :pso_shortcut_straight_max_inserted => :refinement_straight_max_inserted,
    :pso_shortcut_straight_clearance_margin => :refinement_straight_clearance_margin_m,
)

function _rpo_pso_normalize_kwargs(kwargs)
    pairs = Pair{Symbol, Any}[]
    for (key, value) in kwargs
        push!(pairs, get(RPO_PSO_CONFIG_ALIASES, key, key) => value)
    end
    return (; pairs...)
end

function RPOPSOConfig(configurator::RPOPSOConfigurator; kwargs...)
    swarm = configurator.swarm
    objective = configurator.objective
    adaptive = configurator.adaptive
    cull = configurator.cull
    schedule = configurator.schedule
    stagnation = configurator.stagnation
    probe = configurator.probe
    reexplore = configurator.reexplore
    gif = configurator.gif
    refinement = configurator.refinement
    retiming = configurator.retiming
    cfg = RPOPSOConfig(;
        n_waypoints=swarm.n_waypoints,
        n_particles=swarm.n_particles,
        n_iters=swarm.n_iters,
        w_len=objective.w_len,
        w_detour=objective.w_detour,
        w_dir=objective.w_dir,
        w_obs=objective.w_obs,
        w_fuel=objective.w_fuel,
        fuel_proxy_mode=objective.fuel_proxy_mode,
        w_inertia=objective.w_inertia,
        c1=objective.c1,
        c2=objective.c2,
        spread_scale=swarm.spread_scale,
        search_margin_m=swarm.search_margin_m,
        sample_ds_m=swarm.sample_ds_m,
        curve_type=swarm.curve_type,
        detour_ratio_ref=objective.detour_ratio_ref,
        dist_ref_m=objective.dist_ref_m,
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
        cull_enable=cull.enabled,
        cull_fraction_max=cull.fraction_max,
        cull_start_iter=cull.start_iter,
        cull_noise_scale=cull.noise_scale,
        cull_corner_pull=cull.corner_pull,
        cull_obstacle_pull=cull.obstacle_pull,
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
        single_run_gif_enable=gif.single_run_enable,
        gif_top_k=gif.top_k,
        gif_fps=gif.fps,
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
    return isempty(kwargs) ? validate_rpo_pso_config(cfg) : rpo_pso_config(cfg; kwargs...)
end

function rpo_pso_config(base::RPOPSOConfig=RPOPSOConfig(); kwargs...)
    values = merge(_rpo_pso_config_tuple(base), _rpo_pso_normalize_kwargs(kwargs))
    return validate_rpo_pso_config(RPOPSOConfig(; values...))
end

rpo_pso_config(configurator::RPOPSOConfigurator; kwargs...) = RPOPSOConfig(configurator; kwargs...)

function validate_rpo_pso_config(cfg::RPOPSOConfig)
    cfg.n_waypoints >= 0 || throw(ArgumentError("n_waypoints must be nonnegative."))
    cfg.n_particles > 0 || throw(ArgumentError("n_particles must be positive."))
    cfg.n_iters >= 0 || throw(ArgumentError("n_iters must be nonnegative."))
    cfg.sample_ds_m > 0.0 || throw(ArgumentError("sample_ds_m must be positive."))
    cfg.curve_type in (:bezier, :polyline) ||
        throw(ArgumentError("curve_type must be :bezier or :polyline."))
    cfg.fuel_proxy_mode in (:retimed, :nominal) ||
        throw(ArgumentError("fuel_proxy_mode must be :retimed or :nominal."))
    cfg.w_detour >= 0.0 || throw(ArgumentError("w_detour must be nonnegative."))
    cfg.w_dir >= 0.0 || throw(ArgumentError("w_dir must be nonnegative."))
    cfg.detour_ratio_ref > 0.0 || throw(ArgumentError("detour_ratio_ref must be positive."))
    cfg.dist_ref_m > 0.0 || throw(ArgumentError("dist_ref_m must be positive."))
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
    0.0 <= cfg.cull_fraction_max <= 1.0 ||
        throw(ArgumentError("cull_fraction_max must be between 0 and 1."))
    cfg.cull_start_iter >= 0 || throw(ArgumentError("cull_start_iter must be nonnegative."))
    cfg.cull_noise_scale >= 0.0 || throw(ArgumentError("cull_noise_scale must be nonnegative."))
    cfg.cull_corner_pull >= 0.0 || throw(ArgumentError("cull_corner_pull must be nonnegative."))
    cfg.cull_obstacle_pull >= 0.0 || throw(ArgumentError("cull_obstacle_pull must be nonnegative."))
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
    cfg.gif_top_k >= 0 || throw(ArgumentError("gif_top_k must be nonnegative."))
    cfg.gif_fps > 0 || throw(ArgumentError("gif_fps must be positive."))
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
