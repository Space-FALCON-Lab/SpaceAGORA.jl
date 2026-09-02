function save_checkpoint(path::AbstractString, learner::A3CLearner; manifest=nothing)
    mkpath(dirname(path))
    payload = Dict(
        :algorithm => :a3c,
        :actor => cpu_network(learner.actor),
        :critic => cpu_network(learner.critic),
        :actor_optimizer => cpu_adam_state(learner.actor_optimizer),
        :critic_optimizer => cpu_adam_state(learner.critic_optimizer),
        :config => learner.config,
        :action_config => learner.action_config,
        :action_space => learner.action_config.mode,
        :global_step => learner.global_step,
        :train_steps => learner.train_steps,
        :last_loss => learner.last_loss,
        :mean_training_loss => mean_training_loss(learner),
        :training_loss_sum => learner.loss_sum,
        :training_loss_count => learner.loss_count,
        :last_policy_loss => learner.last_policy_loss,
        :last_entropy => learner.last_entropy,
        :last_value_loss => learner.last_value_loss,
        :last_explained_variance => learner.last_explained_variance,
        :last_actor_gradient_norm => learner.last_actor_gradient_norm,
        :last_critic_gradient_norm => learner.last_critic_gradient_norm,
        :policy_version => learner.policy_version,
        :mean_policy_lag => mean_policy_lag(learner),
        :max_observed_policy_lag => learner.max_observed_policy_lag,
        :dropped_stale_updates => learner.dropped_stale_updates,
        :device => training_device_name(learner.device),
        :action_table => copy(PAPER_ACTIONS_MPS),
        :action_bounds_mps => uses_continuous_actions(learner.action_config) ?
                              (CONTINUOUS_ACTION_LOW_MPS, CONTINUOUS_ACTION_HIGH_MPS) :
                              nothing,
        :manifest => _checkpoint_manifest_payload(manifest),
    )
    serialize(path, payload)
    return path
end
