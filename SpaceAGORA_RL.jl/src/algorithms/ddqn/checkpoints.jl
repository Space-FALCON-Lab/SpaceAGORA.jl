function cpu_adam_state(opt::AdamState)
    return AdamState(cpu_network(opt.m), cpu_network(opt.v), opt.t, opt.beta1, opt.beta2,
                     opt.epsilon, opt.learning_rate)
end

function save_checkpoint(path::AbstractString, learner::DDQNLearner; manifest=nothing)
    mkpath(dirname(path))
    algorithm = manifest isa RunManifest ? manifest.algorithm : :ddqn
    payload = Dict(
        :algorithm => algorithm,
        :online => cpu_network(learner.online),
        :target => cpu_network(learner.target),
        :optimizer => cpu_adam_state(learner.optimizer),
        :config => learner.config,
        :schedule => learner.schedule,
        :global_step => learner.global_step,
        :train_steps => learner.train_steps,
        :last_loss => learner.last_loss,
        :device => training_device_name(learner.device),
        :action_table => copy(PAPER_ACTIONS_MPS),
        :manifest => manifest,
    )
    serialize(path, payload)
    return path
end

function save_checkpoint(path::AbstractString, learner::A2CLearner; manifest=nothing)
    mkpath(dirname(path))
    payload = Dict(
        :algorithm => :a2c,
        :actor => cpu_network(learner.actor),
        :critic => cpu_network(learner.critic),
        :actor_optimizer => cpu_adam_state(learner.actor_optimizer),
        :critic_optimizer => cpu_adam_state(learner.critic_optimizer),
        :config => learner.config,
        :global_step => learner.global_step,
        :train_steps => learner.train_steps,
        :last_loss => learner.last_loss,
        :last_policy_loss => learner.last_policy_loss,
        :last_entropy => learner.last_entropy,
        :last_value_loss => learner.last_value_loss,
        :device => training_device_name(learner.device),
        :action_table => copy(PAPER_ACTIONS_MPS),
        :manifest => manifest,
    )
    serialize(path, payload)
    return path
end

load_checkpoint(path::AbstractString) = deserialize(path)
