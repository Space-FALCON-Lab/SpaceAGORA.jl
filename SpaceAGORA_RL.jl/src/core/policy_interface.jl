abstract type AbstractPolicy end

"""
    policy_action_index(policy, config, state, observation, rng)

Return a 1-based discrete action index, or an `AerobrakingAction` for
continuous baseline policies.
"""
function policy_action_index end

function _resolve_policy_action(selected)
    if selected isa AerobrakingAction
        return selected.index, selected
    end
    action_index = Int(selected)
    return action_index, action_from_index(action_index)
end
