abstract type AbstractMDP end

"""
    observe_state(config, state)

Return the task-specific raw observation for a backend decision state.
"""
function observe_state end

"""
    normalize_observation(obs, bounds)

Return the learner-facing normalized observation vector.
"""
function normalize_observation end
