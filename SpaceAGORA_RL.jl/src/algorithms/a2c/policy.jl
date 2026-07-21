struct GreedyA2CPolicy <: AbstractPolicy
    actor::QNetwork
end

function policy_action_index(policy::GreedyA2CPolicy, config::AerobrakingScenarioConfig,
                             state, observation::PaperObservation, rng::AbstractRNG)
    normalized = normalize_observation(observation, config.normalization_bounds)
    return actor_action(policy.actor, normalized, rng; test=true)
end

function load_trained_a2c_policy(checkpoint_path::AbstractString)
    payload = load_checkpoint(checkpoint_path)
    haskey(payload, :actor) || throw(ArgumentError("checkpoint does not contain an :actor network"))
    payload[:actor] isa QNetwork || throw(ArgumentError("checkpoint :actor entry is not a QNetwork"))
    return GreedyA2CPolicy(payload[:actor])
end

