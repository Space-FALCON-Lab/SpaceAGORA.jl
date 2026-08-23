struct RPOEditorAction
    kind::Symbol
    slot::Int
    axis::Int
    direction::Int
    scale::Int
end

function rpo_editor_actions(config::RPOHyPRRLConfig)
    actions = RPOEditorAction[]
    for slot in 1:config.max_translation_waypoints, axis in 1:3,
        direction in (-1, 1), scale in 1:length(config.translation_step_m)
        push!(actions, RPOEditorAction(:move_translation, slot, axis, direction, scale))
    end
    for segment in 1:(config.max_translation_waypoints + 1)
        push!(actions, RPOEditorAction(:insert_translation, segment, 0, 0, 0))
    end
    for slot in 1:config.max_translation_waypoints
        push!(actions, RPOEditorAction(:delete_translation, slot, 0, 0, 0))
    end
    for slot in 1:config.max_attitude_waypoints, axis in 1:3,
        direction in (-1, 1), scale in 1:length(config.attitude_step_rad)
        push!(actions, RPOEditorAction(:rotate_attitude, slot, axis, direction, scale))
    end
    for slot in 1:config.max_attitude_waypoints, direction in (-1, 1)
        push!(actions, RPOEditorAction(:shift_attitude, slot, 0, direction, 0))
    end
    for segment in 1:(config.max_attitude_waypoints + 1)
        push!(actions, RPOEditorAction(:insert_attitude, segment, 0, 0, 0))
    end
    for slot in 1:config.max_attitude_waypoints
        push!(actions, RPOEditorAction(:delete_attitude, slot, 0, 0, 0))
    end
    push!(actions, RPOEditorAction(:stop, 0, 0, 0, 0))
    return actions
end

action_count(config::RPOHyPRRLConfig) = length(rpo_editor_actions(config))

function valid_action_mask(mdp::RPOHyPRRLMDP, state::RPOHyPRRLState)
    n_translation = size(state.control_points_rtn, 2) - 2
    n_attitude = length(state.attitude_progress) - 2
    n_translation_segments = n_translation + 1
    n_attitude_segments = n_attitude + 1
    mask = BitVector(undef, action_count(mdp.config))
    for (index, action) in pairs(rpo_editor_actions(mdp.config))
        mask[index] = if action.kind == :move_translation
            action.slot <= n_translation
        elseif action.kind == :insert_translation
            n_translation < mdp.config.max_translation_waypoints &&
                action.slot <= n_translation_segments
        elseif action.kind == :delete_translation
            action.slot <= n_translation
        elseif action.kind == :rotate_attitude || action.kind == :shift_attitude
            action.slot <= n_attitude
        elseif action.kind == :insert_attitude
            n_attitude < mdp.config.max_attitude_waypoints &&
                action.slot <= n_attitude_segments
        elseif action.kind == :delete_attitude
            action.slot <= n_attitude
        else
            action.kind == :stop
        end
    end
    return mask
end

function _path_frame(points::AbstractMatrix{<:Real}, internal_slot::Int)
    column = internal_slot + 1
    tangent = Vector{Float64}(points[:, min(column + 1, size(points, 2))] -
                              points[:, max(column - 1, 1)])
    tangent ./= max(norm(tangent), eps(Float64))
    reference = abs(tangent[3]) < 0.9 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    normal = cross(reference, tangent)
    normal ./= max(norm(normal), eps(Float64))
    binormal = cross(tangent, normal)
    return (tangent, normal, binormal)
end
