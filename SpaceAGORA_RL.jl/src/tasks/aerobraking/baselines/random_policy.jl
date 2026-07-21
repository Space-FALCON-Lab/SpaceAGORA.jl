struct RandomActionPolicy <: AbstractPolicy end

function policy_action_index(::RandomActionPolicy, config, state, observation, rng::AbstractRNG)
    return rand(rng, 1:action_count())
end
