function save_checkpoint(path::AbstractString, learner::TD3Learner; manifest=nothing)
    mkpath(dirname(path))
    payload = Dict(
        :algorithm => :td3,
        :actor => cpu_network(learner.actor),
        :target_actor => cpu_network(learner.target_actor),
        :critic1 => cpu_network(learner.critic1),
        :critic2 => cpu_network(learner.critic2),
        :target_critic1 => cpu_network(learner.target_critic1),
        :target_critic2 => cpu_network(learner.target_critic2),
        :actor_optimizer => cpu_adam_state(learner.actor_optimizer),
        :critic1_optimizer => cpu_adam_state(learner.critic1_optimizer),
        :critic2_optimizer => cpu_adam_state(learner.critic2_optimizer),
        :config => learner.config,
        :global_step => learner.global_step,
        :train_steps => learner.train_steps,
        :actor_updates => learner.actor_updates,
        :last_actor_loss => learner.last_actor_loss,
        :last_critic1_loss => learner.last_critic1_loss,
        :last_critic2_loss => learner.last_critic2_loss,
        :mean_training_loss => mean_training_loss(learner),
        :last_actor_gradient_norm => learner.last_actor_gradient_norm,
        :last_critic1_gradient_norm => learner.last_critic1_gradient_norm,
        :last_critic2_gradient_norm => learner.last_critic2_gradient_norm,
        :device => training_device_name(learner.device),
        :action_bounds_mps => (CONTINUOUS_ACTION_LOW_MPS, CONTINUOUS_ACTION_HIGH_MPS),
        :manifest => _checkpoint_manifest_payload(manifest),
    )
    serialize(path, payload)
    return path
end
