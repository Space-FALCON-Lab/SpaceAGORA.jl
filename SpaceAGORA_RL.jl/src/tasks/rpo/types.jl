Base.@kwdef struct RPOHyPRRLConfig
    # Internal Bezier control waypoints; the start and goal are fixed endpoints.
    max_translation_waypoints::Int = 20
    max_attitude_waypoints::Int = 6
    max_edits::Int = 64
    translation_step_m::NTuple{2, Float64} = (0.10, 0.40)
    attitude_step_rad::NTuple{2, Float64} = (deg2rad(2.0), deg2rad(10.0))
    attitude_time_step::Float64 = 0.05
    safe_distance_m::Float64 = 0.5
    coordinate_scale_m::Float64 = 20.0
    reward_improvement_scale::Float64 = 1.0
    invalid_action_penalty::Float64 = 0.25
    infeasible_edit_penalty::Float64 = 1.0
    edit_penalty::Float64 = 0.01
    completion_bonus::Float64 = 0.25
    fuel_weight::Float64 = 1.0
    duration_weight::Float64 = 0.0
    allocation_error_weight::Float64 = 0.0
    wheel_weight::Float64 = 0.01
    mass_kg::Float64 = 5.2
    g0_mps2::Float64 = 9.80665
    thruster_directions_body::Matrix{Float64} = [
         1.0 -1.0  0.0  0.0  0.0  0.0
         0.0  0.0  1.0 -1.0  0.0  0.0
         0.0  0.0  0.0  0.0  1.0 -1.0
    ]
    thruster_max_thrust_n::Vector{Float64} = fill(0.05, 6)
    thruster_isp_s::Vector{Float64} = fill(60.0, 6)
    thruster_min_firing_time_s::Float64 = 0.010
    reaction_wheel_inertia_kgm2::Vector{Float64} = [0.035, 0.035, 0.025]
    reaction_wheel_kp::Float64 = 0.020
    reaction_wheel_kd::Float64 = 0.050
    reaction_wheel_max_torque_nm::Float64 = 0.002
    reaction_wheel_max_momentum_nms::Float64 = 0.020
end

function _validate_rpo_hypr_rl_config(config::RPOHyPRRLConfig)
    config.max_translation_waypoints >= 0 ||
        throw(ArgumentError("max_translation_waypoints must be nonnegative"))
    config.max_attitude_waypoints >= 0 ||
        throw(ArgumentError("max_attitude_waypoints must be nonnegative"))
    config.max_edits > 0 || throw(ArgumentError("max_edits must be positive"))
    config.reward_improvement_scale > 0.0 ||
        throw(ArgumentError("reward_improvement_scale must be positive"))
    config.invalid_action_penalty >= 0.0 ||
        throw(ArgumentError("invalid_action_penalty must be nonnegative"))
    config.infeasible_edit_penalty > 0.0 ||
        throw(ArgumentError("infeasible_edit_penalty must be positive"))
    config.edit_penalty >= 0.0 ||
        throw(ArgumentError("edit_penalty must be nonnegative"))
    config.completion_bonus >= 0.0 ||
        throw(ArgumentError("completion_bonus must be nonnegative"))
    all(>=(0.0), (
        config.fuel_weight, config.duration_weight,
        config.allocation_error_weight, config.wheel_weight,
    )) || throw(ArgumentError("objective weights must be nonnegative"))
    all(>(0.0), config.translation_step_m) ||
        throw(ArgumentError("translation steps must be positive"))
    all(>(0.0), config.attitude_step_rad) ||
        throw(ArgumentError("attitude steps must be positive"))
    config.safe_distance_m >= 0.0 ||
        throw(ArgumentError("safe_distance_m must be nonnegative"))
    config.mass_kg > 0.0 || throw(ArgumentError("mass_kg must be positive"))
    config.g0_mps2 > 0.0 || throw(ArgumentError("g0_mps2 must be positive"))
    size(config.thruster_directions_body) == (3, 6) ||
        throw(DimensionMismatch("thruster_directions_body must be 3x6"))
    length(config.thruster_max_thrust_n) == 6 ||
        throw(DimensionMismatch("six maximum-thrust values are required"))
    length(config.thruster_isp_s) == 6 ||
        throw(DimensionMismatch("six specific-impulse values are required"))
    all(>(0.0), config.thruster_isp_s) ||
        throw(ArgumentError("thruster specific impulse must be positive"))
    all(>=(0.0), config.thruster_max_thrust_n) ||
        throw(ArgumentError("maximum thrust must be nonnegative"))
    config.thruster_min_firing_time_s >= 0.0 ||
        throw(ArgumentError("minimum firing time must be nonnegative"))
    length(config.reaction_wheel_inertia_kgm2) == 3 ||
        throw(DimensionMismatch("three reaction-wheel inertia values are required"))
    all(>(0.0), config.reaction_wheel_inertia_kgm2) ||
        throw(ArgumentError("reaction-wheel inertia must be positive"))
    config.reaction_wheel_max_torque_nm >= 0.0 ||
        throw(ArgumentError("reaction-wheel torque limit must be nonnegative"))
    config.reaction_wheel_max_momentum_nms > 0.0 ||
        throw(ArgumentError("reaction-wheel momentum limit must be positive"))
    return config
end

Base.@kwdef struct RPOHyPRRLScenario
    start_rtn::Vector{Float64}
    goal_rtn::Vector{Float64}
    geometry::Any
    pso_config::Any = nothing
    tracking_settings::Any = nothing
    rrt_settings::Any = nothing
    seed_path_rtn::Union{Nothing, Matrix{Float64}} = nothing
    initial_attitude_rtn_to_body::Vector{Float64} = [0.0, 0.0, 0.0, 1.0]
    final_attitude_rtn_to_body::Vector{Float64} = [0.0, 0.0, 0.0, 1.0]
end

struct RPOHyPRRLEvaluation
    feasible::Bool
    objective::Float64
    path_cost::Float64
    propellant_used_kg::Float64
    duration_s::Float64
    allocation_error_impulse_mps::Float64
    thruster_saturation_fraction::Float64
    wheel_energy_j::Float64
    wheel_peak_momentum_nms::Float64
    min_clearance_m::Float64
    final_position_error_m::Float64
    t_ref_s::Vector{Float64}
    r_ref_rtn::Matrix{Float64}
    v_ref_rtn::Matrix{Float64}
    q_ref_rtn_to_body::Matrix{Float64}
    thruster_impulse_ns::Vector{Float64}
    diagnostics::NamedTuple
end

function _failed_rpo_evaluation(; reason::Symbol=:evaluation_failed,
                                error_message::AbstractString="")
    return RPOHyPRRLEvaluation(
        false, Inf, Inf, Inf, Inf, Inf, 1.0, Inf, Inf, -Inf, Inf,
        Float64[], zeros(3, 0), zeros(3, 0), zeros(4, 0), zeros(6),
        (reason=reason, error_message=String(error_message)),
    )
end

struct RPOHyPRRLState
    control_points_rtn::Matrix{Float64}
    attitude_progress::Vector{Float64}
    attitude_quaternions::Matrix{Float64}
    evaluation::RPOHyPRRLEvaluation
    best_control_points_rtn::Matrix{Float64}
    best_attitude_progress::Vector{Float64}
    best_attitude_quaternions::Matrix{Float64}
    best_evaluation::RPOHyPRRLEvaluation
    seed_evaluation::RPOHyPRRLEvaluation
    edit_count::Int
    stopped::Bool
end

struct RPOHyPRRLStepResult
    state::RPOHyPRRLState
    reward::Float64
    terminated::Bool
    truncated::Bool
    accepted::Bool
    action_index::Int
    info::NamedTuple
end

struct RPOHyPRRLPlan
    valid::Bool
    path_rtn::Matrix{Float64}
    t_ref_s::Vector{Float64}
    r_ref_rtn::Matrix{Float64}
    v_ref_rtn::Matrix{Float64}
    q_ref_rtn_to_body::Matrix{Float64}
    propellant_used_kg::Float64
    cost::Float64
    seed_cost::Float64
    edits::Int
    stopped::Bool
    diagnostics::NamedTuple
end

struct RPOHyPRRLMDP{F} <: AbstractMDP
    config::RPOHyPRRLConfig
    scenario::RPOHyPRRLScenario
    evaluator::F
end

function RPOHyPRRLMDP(config::RPOHyPRRLConfig, scenario::RPOHyPRRLScenario;
                      evaluator=evaluate_rpo_candidate)
    _validate_rpo_hypr_rl_config(config)
    length(scenario.start_rtn) == 3 || throw(DimensionMismatch("start_rtn must have length 3"))
    length(scenario.goal_rtn) == 3 || throw(DimensionMismatch("goal_rtn must have length 3"))
    return RPOHyPRRLMDP{typeof(evaluator)}(config, scenario, evaluator)
end
