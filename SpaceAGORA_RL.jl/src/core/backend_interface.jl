abstract type AbstractRLBackend end

"""
    reset_scenario(config, rng) -> state

Create a backend-owned decision state for a new episode.
"""
function reset_scenario end

"""
    step_scenario(config, state, action, rng) -> step_result

Advance one RL decision step and return transition data.
"""
function step_scenario end
