struct GreedyDDQNPolicy <: AbstractPolicy
    network::QNetwork
end

const GreedyPRDRLPolicy = GreedyDDQNPolicy

function policy_action_index(policy::GreedyDDQNPolicy, config::AerobrakingScenarioConfig,
                             state, observation::PaperObservation, rng::AbstractRNG)
    normalized = normalize_observation(observation, config.normalization_bounds)
    return argmax(predict_q(policy.network, normalized))
end

function load_trained_ddqn_policy(checkpoint_path::AbstractString)
    payload = load_checkpoint(checkpoint_path)
    haskey(payload, :online) || throw(ArgumentError("checkpoint does not contain an :online Q-network"))
    payload[:online] isa QNetwork || throw(ArgumentError("checkpoint :online entry is not a QNetwork"))
    return GreedyDDQNPolicy(payload[:online])
end

load_trained_pr_drl_policy(checkpoint_path::AbstractString) = load_trained_ddqn_policy(checkpoint_path)
