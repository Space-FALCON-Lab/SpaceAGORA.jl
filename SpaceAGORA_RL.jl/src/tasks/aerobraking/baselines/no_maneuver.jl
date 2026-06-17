struct NoManeuverPolicy <: AbstractPolicy end

function policy_action_index(::NoManeuverPolicy, config, state, observation, rng::AbstractRNG)
    return zero_action_index()
end
