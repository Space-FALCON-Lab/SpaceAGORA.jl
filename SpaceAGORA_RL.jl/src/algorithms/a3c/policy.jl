struct GreedyA3CPolicy <: AbstractPolicy
    actor::QNetwork
    action_config::ActorCriticActionConfig
end

GreedyA3CPolicy(actor::QNetwork) = GreedyA3CPolicy(actor, ActorCriticActionConfig())

function policy_action_index(policy::GreedyA3CPolicy, config::AerobrakingScenarioConfig,
                             state, observation::PaperObservation, rng::AbstractRNG)
    normalized = normalize_observation(observation, config.normalization_bounds)
    return actor_critic_action(
        policy.actor,
        policy.action_config,
        normalized,
        rng;
        test=true,
    )
end

function load_trained_a3c_policy(checkpoint_path::AbstractString)
    payload = load_checkpoint(checkpoint_path)
    get(payload, :algorithm, nothing) == :a3c ||
        throw(ArgumentError("checkpoint is not an A3C checkpoint"))
    haskey(payload, :actor) || throw(ArgumentError("checkpoint does not contain an :actor network"))
    payload[:actor] isa QNetwork || throw(ArgumentError("checkpoint :actor entry is not a QNetwork"))
    action_config = get(payload, :action_config, ActorCriticActionConfig())
    action_config isa ActorCriticActionConfig ||
        throw(ArgumentError("checkpoint :action_config entry is invalid"))
    return GreedyA3CPolicy(payload[:actor], action_config)
end
