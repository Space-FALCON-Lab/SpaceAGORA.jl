struct GreedyTD3Policy <: AbstractPolicy
    actor::QNetwork
end

function policy_action_index(policy::GreedyTD3Policy, config::AerobrakingScenarioConfig,
                             state, observation::PaperObservation, rng::AbstractRNG)
    normalized_observation = normalize_observation(observation, config.normalization_bounds)
    normalized_action = only(td3_actor_output(policy.actor, normalized_observation))
    return _continuous_action_from_normalized(normalized_action)
end

function load_trained_td3_policy(checkpoint_path::AbstractString)
    payload = load_checkpoint(checkpoint_path)
    get(payload, :algorithm, nothing) == :td3 ||
        throw(ArgumentError("checkpoint is not a TD3 checkpoint"))
    haskey(payload, :actor) || throw(ArgumentError("TD3 checkpoint does not contain an :actor network"))
    payload[:actor] isa QNetwork || throw(ArgumentError("TD3 checkpoint :actor entry is not a QNetwork"))
    return GreedyTD3Policy(payload[:actor])
end
