abstract type AbstractPolicy end

"""
    policy_action_index(policy, config, state, observation, rng)

Return a 1-based discrete action index.
"""
function policy_action_index end
