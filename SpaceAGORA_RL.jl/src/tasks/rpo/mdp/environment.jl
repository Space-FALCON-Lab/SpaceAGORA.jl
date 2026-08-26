function _initial_attitude_knots(scenario::RPOHyPRRLScenario)
    progress = [0.0, 1.0]
    quaternions = hcat(
        _quat_normalize(scenario.initial_attitude_rtn_to_body),
        _quat_normalize(scenario.final_attitude_rtn_to_body),
    )
    return progress, quaternions
end

function reset_scenario(mdp::RPOHyPRRLMDP, rng::AbstractRNG=Random.default_rng())
    seed_path = rpo_hypr_rl_seed_path(mdp.scenario, mdp.config; rng=rng)
    progress, quaternions = _initial_attitude_knots(mdp.scenario)
    evaluation = mdp.evaluator(
        mdp.scenario, mdp.config, seed_path, progress, quaternions,
    )
    return RPOHyPRRLState(
        seed_path, progress, quaternions, evaluation,
        copy(seed_path), copy(progress), copy(quaternions), evaluation,
        evaluation, 0, false,
    )
end

function _insert_matrix_column(matrix::AbstractMatrix, after_column::Int, value)
    return hcat(matrix[:, 1:after_column], value, matrix[:, (after_column + 1):end])
end

function _delete_matrix_column(matrix::AbstractMatrix, column::Int)
    keep = [index for index in axes(matrix, 2) if index != column]
    return Matrix{Float64}(matrix[:, keep])
end

function _rpo_infeasibility_score(evaluation::RPOHyPRRLEvaluation,
                                  config::RPOHyPRRLConfig)
    evaluation.feasible && return 0.0
    components = get(evaluation.diagnostics, :path_components, nothing)
    components === nothing && return Inf
    penalty = Float64(get(
        evaluation.diagnostics,
        :clearance_penalty,
        get(components, :total, Inf),
    ))
    return isfinite(penalty) ? max(0.0, penalty) : Inf
end

function _rpo_infeasibility_improvement(before::RPOHyPRRLEvaluation,
                                        after::RPOHyPRRLEvaluation,
                                        config::RPOHyPRRLConfig)
    before_score = _rpo_infeasibility_score(before, config)
    after_score = _rpo_infeasibility_score(after, config)
    limit = config.infeasible_edit_penalty
    if isfinite(before_score) && isfinite(after_score)
        return clamp(before_score - after_score, -limit, limit)
    elseif !isfinite(before_score) && isfinite(after_score)
        return limit
    elseif isfinite(before_score) && !isfinite(after_score)
        return -limit
    end
    return 0.0
end

function _rpo_bounded_objective_improvement(before::Real, after::Real,
                                            config::RPOHyPRRLConfig;
                                            allow_negative::Bool=true)
    if !isfinite(before) || !isfinite(after)
        return 0.0
    end
    scaled = (Float64(before) - Float64(after)) /
        max(config.reward_improvement_scale, 1.0e-9)
    lower = allow_negative ? -config.infeasible_edit_penalty : 0.0
    return clamp(scaled, lower, config.infeasible_edit_penalty)
end

function _rpo_terminal_reward(seed::RPOHyPRRLEvaluation,
                              best::RPOHyPRRLEvaluation,
                              config::RPOHyPRRLConfig)
    best.feasible || return -config.infeasible_edit_penalty
    seed.feasible || return config.infeasible_edit_penalty
    improvement = _rpo_bounded_objective_improvement(
        seed.objective, best.objective, config; allow_negative=false,
    )
    return config.completion_bonus * improvement
end

function _apply_rpo_editor_action(state::RPOHyPRRLState,
                                  action::RPOEditorAction,
                                  config::RPOHyPRRLConfig)
    points = copy(state.control_points_rtn)
    progress = copy(state.attitude_progress)
    quaternions = copy(state.attitude_quaternions)
    if action.kind == :move_translation
        frame = _path_frame(points, action.slot)
        delta = action.direction * config.translation_step_m[action.scale] .* frame[action.axis]
        points[:, action.slot + 1] .+= delta
    elseif action.kind == :insert_translation
        midpoint = 0.5 .* (points[:, action.slot] .+ points[:, action.slot + 1])
        points = _insert_matrix_column(points, action.slot, midpoint)
    elseif action.kind == :delete_translation
        points = _delete_matrix_column(points, action.slot + 1)
    elseif action.kind == :rotate_attitude
        column = action.slot + 1
        rotation = _axis_angle_quaternion(
            action.axis,
            action.direction * config.attitude_step_rad[action.scale],
        )
        quaternions[:, column] .= _quat_normalize(
            _quat_multiply(rotation, view(quaternions, :, column)),
        )
    elseif action.kind == :shift_attitude
        index = action.slot + 1
        lower = progress[index - 1] + 1.0e-3
        upper = progress[index + 1] - 1.0e-3
        progress[index] = clamp(
            progress[index] + action.direction * config.attitude_time_step,
            lower, upper,
        )
    elseif action.kind == :insert_attitude
        lower = action.slot
        upper = lower + 1
        new_progress = 0.5 * (progress[lower] + progress[upper])
        new_quaternion = _quat_slerp(
            view(quaternions, :, lower), view(quaternions, :, upper), 0.5,
        )
        insert!(progress, upper, new_progress)
        quaternions = _insert_matrix_column(quaternions, lower, new_quaternion)
    elseif action.kind == :delete_attitude
        deleteat!(progress, action.slot + 1)
        quaternions = _delete_matrix_column(quaternions, action.slot + 1)
    end
    return points, progress, quaternions
end

function step_scenario(mdp::RPOHyPRRLMDP, state::RPOHyPRRLState,
                       action_index::Integer,
                       rng::AbstractRNG=Random.default_rng())
    state.stopped && throw(ArgumentError("cannot edit a stopped HyPR-RL state"))
    actions = rpo_editor_actions(mdp.config)
    index = Int(action_index)
    1 <= index <= length(actions) || throw(BoundsError(actions, index))
    action = actions[index]
    mask = valid_action_mask(mdp, state)
    next_edit_count = state.edit_count + 1
    if !mask[index]
        next_state = RPOHyPRRLState(
            state.control_points_rtn, state.attitude_progress,
            state.attitude_quaternions, state.evaluation,
            state.best_control_points_rtn, state.best_attitude_progress,
            state.best_attitude_quaternions, state.best_evaluation,
            state.seed_evaluation, next_edit_count, false,
        )
        truncated = next_edit_count >= mdp.config.max_edits
        return RPOHyPRRLStepResult(
            next_state, -mdp.config.invalid_action_penalty, false, truncated,
            false, index, (action=action, reason=:masked_action),
        )
    end
    if action.kind == :stop
        reward = _rpo_terminal_reward(
            state.seed_evaluation, state.best_evaluation, mdp.config,
        )
        next_state = RPOHyPRRLState(
            state.control_points_rtn, state.attitude_progress,
            state.attitude_quaternions, state.evaluation,
            state.best_control_points_rtn, state.best_attitude_progress,
            state.best_attitude_quaternions, state.best_evaluation,
            state.seed_evaluation, next_edit_count, true,
        )
        return RPOHyPRRLStepResult(
            next_state, reward, true, false, true, index,
            (action=action, reason=:stopped),
        )
    end

    points, progress, quaternions = _apply_rpo_editor_action(state, action, mdp.config)
    candidate = mdp.evaluator(mdp.scenario, mdp.config, points, progress, quaternions)
    if !isfinite(candidate.objective)
        next_state = RPOHyPRRLState(
            state.control_points_rtn, state.attitude_progress,
            state.attitude_quaternions, state.evaluation,
            state.best_control_points_rtn, state.best_attitude_progress,
            state.best_attitude_quaternions, state.best_evaluation,
            state.seed_evaluation, next_edit_count, false,
        )
        truncated = next_edit_count >= mdp.config.max_edits
        return RPOHyPRRLStepResult(
            next_state, -mdp.config.infeasible_edit_penalty, false, truncated,
            false, index, (action=action, reason=:evaluation_failed),
        )
    end

    improvement = if state.evaluation.feasible && candidate.feasible
        state.evaluation.objective - candidate.objective
    elseif !state.evaluation.feasible && candidate.feasible
        mdp.config.infeasible_edit_penalty
    else
        _rpo_infeasibility_improvement(
            state.evaluation, candidate, mdp.config,
        )
    end
    is_best = if candidate.feasible
        !state.best_evaluation.feasible ||
            candidate.objective + 1.0e-12 < state.best_evaluation.objective
    elseif !state.best_evaluation.feasible
        _rpo_infeasibility_score(candidate, mdp.config) + 1.0e-12 <
            _rpo_infeasibility_score(state.best_evaluation, mdp.config)
    else
        false
    end
    best_points = is_best ? copy(points) : state.best_control_points_rtn
    best_progress = is_best ? copy(progress) : state.best_attitude_progress
    best_quaternions = is_best ? copy(quaternions) : state.best_attitude_quaternions
    best_evaluation = is_best ? candidate : state.best_evaluation
    next_state = RPOHyPRRLState(
        points, progress, quaternions, candidate,
        best_points, best_progress, best_quaternions, best_evaluation,
        state.seed_evaluation, next_edit_count, false,
    )
    bounded_improvement = clamp(
        improvement / max(mdp.config.reward_improvement_scale, 1.0e-9),
        -mdp.config.infeasible_edit_penalty,
        mdp.config.infeasible_edit_penalty,
    )
    reward = bounded_improvement - mdp.config.edit_penalty
    truncated = next_edit_count >= mdp.config.max_edits
    return RPOHyPRRLStepResult(
        next_state, reward, false, truncated, true, index,
        (
            action=action,
            reason=candidate.feasible ? :accepted : :accepted_infeasible,
            new_best=is_best,
        ),
    )
end

function rpo_hypr_rl_observation_dim(config::RPOHyPRRLConfig)
    return 6 + 8 * config.max_translation_waypoints +
           6 * config.max_attitude_waypoints + 9 + 3
end

function _rpo_point_geometry_features(point, geometry)
    geometry === nothing && return (clearance=0.0, normal=zeros(3))
    try
        modules = _spaceagora_rpo_modules()
        clearance = Float64(Base.invokelatest(
            getproperty(modules.navigation, :rpo_clearance_distance_to_station),
            point, geometry,
        ))
        normal = Vector{Float64}(Base.invokelatest(
            getproperty(modules.navigation, :rpo_surface_normal_from_pointcloud),
            point, geometry,
        ))
        return (clearance=clearance, normal=normal)
    catch
        return (clearance=0.0, normal=zeros(3))
    end
end

function observe_state(mdp::RPOHyPRRLMDP, state::RPOHyPRRLState)
    config = mdp.config
    scale = max(config.coordinate_scale_m, 1.0e-9)
    values = Float32[]
    append!(values, Float32.(mdp.scenario.start_rtn ./ scale))
    append!(values, Float32.(mdp.scenario.goal_rtn ./ scale))
    n_translation = size(state.control_points_rtn, 2) - 2
    for slot in 1:config.max_translation_waypoints
        if slot <= n_translation
            point = view(state.control_points_rtn, :, slot + 1)
            geometry_features = _rpo_point_geometry_features(
                point, mdp.scenario.geometry,
            )
            append!(values, Float32.(point ./ scale))
            push!(values, Float32(geometry_features.clearance / scale))
            append!(values, Float32.(geometry_features.normal))
            push!(values, 1.0f0)
        else
            append!(values, zeros(Float32, 8))
        end
    end
    n_attitude = length(state.attitude_progress) - 2
    for slot in 1:config.max_attitude_waypoints
        if slot <= n_attitude
            index = slot + 1
            push!(values, Float32(state.attitude_progress[index]))
            append!(values, Float32.(view(state.attitude_quaternions, :, index)))
            push!(values, 1.0f0)
        else
            append!(values, zeros(Float32, 6))
        end
    end
    evaluation = state.evaluation
    metrics = (
        log1p(min(evaluation.objective, 1.0e6)) / log(1.0e6 + 1.0),
        log1p(min(evaluation.path_cost, 1.0e6)) / log(1.0e6 + 1.0),
        min(evaluation.propellant_used_kg / 0.2, 10.0),
        min(evaluation.duration_s / 300.0, 10.0),
        min(evaluation.allocation_error_impulse_mps, 10.0),
        evaluation.thruster_saturation_fraction,
        min(evaluation.wheel_peak_momentum_nms /
            max(config.reaction_wheel_max_momentum_nms, 1.0e-9), 10.0),
        clamp(evaluation.min_clearance_m / scale, -10.0, 10.0),
        min(evaluation.final_position_error_m / scale, 10.0),
    )
    append!(values, Float32.(metrics))
    push!(values, Float32(n_translation / max(config.max_translation_waypoints, 1)))
    push!(values, Float32(n_attitude / max(config.max_attitude_waypoints, 1)))
    push!(values, Float32(state.edit_count / config.max_edits))
    length(values) == rpo_hypr_rl_observation_dim(config) ||
        error("internal HyPR-RL observation dimension mismatch")
    return values
end

normalize_observation(observation::AbstractVector{<:Real}, ::RPOHyPRRLConfig) =
    Float32.(observation)
