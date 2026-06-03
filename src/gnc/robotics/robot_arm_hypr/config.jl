"""Spherical workspace obstacle used by robot-arm HYPR clearance checks."""
struct RobotArmSphereObstacle
    center::SVector{3, Float64}
    radius_m::Float64
end

"""Spherical workspace obstacle used by robot-arm HYPR clearance checks."""
RobotArmSphereObstacle(center, radius_m::Real) =
    RobotArmSphereObstacle(SVector{3, Float64}(center), Float64(radius_m))

"""Configuration for the robot-arm HYPR planner, cost function, and retiming."""
Base.@kwdef struct RobotArmHYPRConfig
    n_waypoints::Int = 4
    n_particles::Int = 96
    n_iters::Int = 55
    n_samples::Int = 72
    curve_type::Symbol = :bezier
    safe_distance_m::Float64 = 0.006
    w_len::Float64 = 1.0
    w_smooth::Float64 = 0.25
    w_obs::Float64 = 600.0
    w_inertia::Float64 = 0.72
    c1::Float64 = 1.35
    c2::Float64 = 1.55
    spread_scale::Float64 = 0.22
    max_velocity_fraction::Float64 = 0.18
    schedule_enable::Bool = true
    schedule_transition_fraction::Float64 = 0.45
    schedule_w_min::Float64 = 0.28
    schedule_w_end_fraction::Float64 = 0.65
    schedule_c1_end_fraction::Float64 = 0.75
    schedule_c2_end_fraction::Float64 = 1.25
    schedule_c_min::Float64 = 0.5
    schedule_c_max::Float64 = 2.5
    cull_enable::Bool = true
    cull_fraction_max::Float64 = 0.30
    cull_start_iter::Int = 8
    cull_noise_scale::Float64 = 0.22
    rrt_warmstart_enable::Bool = false
    rrt_warmstart_iters::Int = 250
    rrt_warmstart_step_size_rad::Float64 = 0.25
    rrt_warmstart_goal_sample_rate::Float64 = 0.05
    rrt_warmstart_collision_sample_ds_rad::Float64 = 0.05
    rrt_warmstart_connect_max_steps::Int = 10_000
    rrt_warmstart_shortcut_iters::Int = 40
    rrt_warmstart_runtime_limit_s::Float64 = Inf
    refinement_enable::Bool = true
    refinement_rounds::Int = 8
    refinement_step_fraction::Float64 = 0.12
    refinement_shrink::Float64 = 0.55
    refinement_min_abs_cost_improvement::Float64 = 1.0e-8
    refinement_min_rel_cost_improvement::Float64 = 1.0e-5
    early_stopping_enable::Bool = true
    early_stopping_min_iters::Int = 16
    early_stopping_patience::Int = 12
    early_stopping_min_abs_improvement::Float64 = 1.0e-8
    early_stopping_min_rel_improvement::Float64 = 1.0e-4
    early_stopping_require_feasible::Bool = false
    retime_enable::Bool = false
    retime_max_joint_velocity_rad_s::Float64 = Inf
    retime_max_joint_acceleration_rad_s2::Float64 = Inf
    retime_reaction_time_s::Float64 = 0.0
    retime_reaction_gain::Float64 = 1.0
    retime_min_scale::Float64 = 1.0
    retime_max_scale::Float64 = 100.0
    retime_max_base_force_n::Float64 = Inf
    retime_max_base_torque_nm::Float64 = Inf
    retime_base_wrench_model_gain::Float64 = 1.0
    retime_base_wrench_margin::Float64 = 1.10
    retime_base_wrench_iters::Int = 6
    retime_cloth_physics_enable::Bool = true
    retime_cloth_dt_s::Float64 = 0.025
    retime_cloth_k_translation_n_m::Float64 = 5.0e3
    retime_cloth_c_translation_n_s_m::Float64 = 30.0
    retime_cloth_k_rotation_n_m_rad::Float64 = 15.0
    retime_cloth_c_rotation_n_m_s_rad::Float64 = 0.5
end

"""Robot-arm HYPR output containing the plan, path, costs, and diagnostics."""
struct RobotArmHYPRResult
    plan::RobotArmPlan
    control_points::Matrix{Float64}
    sampled_path::Matrix{Float64}
    cost::Float64
    components::NamedTuple
    cost_history::Vector{Float64}
    config::RobotArmHYPRConfig
    early_stopped::Bool
    early_stop_iter::Int
    cull_replacements::Int
end

"""Validate robot-arm HYPR configuration dimensions and ranges."""
function _validate_robot_arm_hypr_config(cfg::RobotArmHYPRConfig)
    cfg.n_waypoints >= 0 || throw(ArgumentError("n_waypoints must be nonnegative."))
    cfg.n_particles > 0 || throw(ArgumentError("n_particles must be positive."))
    cfg.n_iters > 0 || throw(ArgumentError("n_iters must be positive."))
    cfg.n_samples >= 2 || throw(ArgumentError("n_samples must be at least 2."))
    cfg.curve_type in (:bezier, :polyline) ||
        throw(ArgumentError("curve_type must be :bezier or :polyline."))
    cfg.safe_distance_m >= 0.0 || throw(ArgumentError("safe_distance_m must be nonnegative."))
    cfg.rrt_warmstart_iters >= 0 ||
        throw(ArgumentError("rrt_warmstart_iters must be nonnegative."))
    cfg.rrt_warmstart_step_size_rad > 0.0 ||
        throw(ArgumentError("rrt_warmstart_step_size_rad must be positive."))
    0.0 <= cfg.rrt_warmstart_goal_sample_rate <= 1.0 ||
        throw(ArgumentError("rrt_warmstart_goal_sample_rate must be between 0 and 1."))
    cfg.rrt_warmstart_collision_sample_ds_rad > 0.0 ||
        throw(ArgumentError("rrt_warmstart_collision_sample_ds_rad must be positive."))
    cfg.rrt_warmstart_connect_max_steps > 0 ||
        throw(ArgumentError("rrt_warmstart_connect_max_steps must be positive."))
    cfg.rrt_warmstart_shortcut_iters >= 0 ||
        throw(ArgumentError("rrt_warmstart_shortcut_iters must be nonnegative."))
    cfg.rrt_warmstart_runtime_limit_s >= 0.0 ||
        throw(ArgumentError("rrt_warmstart_runtime_limit_s must be nonnegative."))
    cfg.retime_max_joint_velocity_rad_s > 0.0 ||
        throw(ArgumentError("retime_max_joint_velocity_rad_s must be positive."))
    cfg.retime_max_joint_acceleration_rad_s2 > 0.0 ||
        throw(ArgumentError("retime_max_joint_acceleration_rad_s2 must be positive."))
    cfg.retime_reaction_time_s >= 0.0 ||
        throw(ArgumentError("retime_reaction_time_s must be nonnegative."))
    cfg.retime_reaction_gain >= 0.0 ||
        throw(ArgumentError("retime_reaction_gain must be nonnegative."))
    cfg.retime_min_scale >= 1.0 ||
        throw(ArgumentError("retime_min_scale must be at least 1.0."))
    cfg.retime_max_scale >= cfg.retime_min_scale ||
        throw(ArgumentError("retime_max_scale must be at least retime_min_scale."))
    cfg.retime_max_base_force_n > 0.0 ||
        throw(ArgumentError("retime_max_base_force_n must be positive."))
    cfg.retime_max_base_torque_nm > 0.0 ||
        throw(ArgumentError("retime_max_base_torque_nm must be positive."))
    cfg.retime_base_wrench_model_gain > 0.0 ||
        throw(ArgumentError("retime_base_wrench_model_gain must be positive."))
    cfg.retime_base_wrench_margin >= 1.0 ||
        throw(ArgumentError("retime_base_wrench_margin must be at least 1.0."))
    cfg.retime_base_wrench_iters >= 0 ||
        throw(ArgumentError("retime_base_wrench_iters must be nonnegative."))
    cfg.refinement_rounds >= 0 ||
        throw(ArgumentError("refinement_rounds must be nonnegative."))
    cfg.refinement_step_fraction > 0.0 ||
        throw(ArgumentError("refinement_step_fraction must be positive."))
    0.0 < cfg.refinement_shrink < 1.0 ||
        throw(ArgumentError("refinement_shrink must be in (0, 1)."))
    cfg.refinement_min_abs_cost_improvement >= 0.0 ||
        throw(ArgumentError("refinement_min_abs_cost_improvement must be nonnegative."))
    cfg.refinement_min_rel_cost_improvement >= 0.0 ||
        throw(ArgumentError("refinement_min_rel_cost_improvement must be nonnegative."))
    cfg.retime_cloth_dt_s > 0.0 ||
        throw(ArgumentError("retime_cloth_dt_s must be positive."))
    cfg.retime_cloth_k_translation_n_m >= 0.0 ||
        throw(ArgumentError("retime_cloth_k_translation_n_m must be nonnegative."))
    cfg.retime_cloth_c_translation_n_s_m >= 0.0 ||
        throw(ArgumentError("retime_cloth_c_translation_n_s_m must be nonnegative."))
    cfg.retime_cloth_k_rotation_n_m_rad >= 0.0 ||
        throw(ArgumentError("retime_cloth_k_rotation_n_m_rad must be nonnegative."))
    cfg.retime_cloth_c_rotation_n_m_s_rad >= 0.0 ||
        throw(ArgumentError("retime_cloth_c_rotation_n_m_s_rad must be nonnegative."))
    return cfg
end
