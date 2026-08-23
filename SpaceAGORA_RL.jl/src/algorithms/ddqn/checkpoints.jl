function cpu_adam_state(opt::AdamState)
    return AdamState(cpu_network(opt.m), cpu_network(opt.v), opt.t, opt.beta1, opt.beta2,
                     opt.epsilon, opt.learning_rate)
end

function _deserialize_a2c_config_field(s::Serialization.AbstractSerializer)
    tag = Int32(read(s.io, UInt8)::UInt8)
    tag == Serialization.UNDEFREF_TAG &&
        throw(ArgumentError("serialized A2CConfig contains an undefined field"))
    return Serialization.handle_deserialize(s, tag)
end

function Serialization.deserialize(s::Serialization.AbstractSerializer, ::Type{A2CConfig})
    prefix = ntuple(_ -> _deserialize_a2c_config_field(s), 7)
    if prefix[7] isa Bool
        suffix = ntuple(_ -> _deserialize_a2c_config_field(s), 7)
        return A2CConfig(prefix..., suffix...)
    elseif prefix[7] isa Float64
        suffix = ntuple(_ -> _deserialize_a2c_config_field(s), 6)
        return A2CConfig(prefix[1:6]..., true, prefix[7], suffix...)
    end
    throw(ArgumentError(
        "serialized A2CConfig has an unsupported seventh field of type $(typeof(prefix[7]))",
    ))
end

_checkpoint_manifest_payload(manifest) =
    manifest isa RunManifest ? manifest_dict(manifest) : manifest

function save_checkpoint(path::AbstractString, learner::DDQNLearner;
                         manifest=nothing, task::Symbol=:aerobraking,
                         action_table=copy(PAPER_ACTIONS_MPS), task_metadata=nothing)
    mkpath(dirname(path))
    algorithm = manifest isa RunManifest ? manifest.algorithm : :ddqn
    payload = Dict(
        :algorithm => algorithm,
        :task => task,
        :online => cpu_network(learner.online),
        :target => cpu_network(learner.target),
        :optimizer => cpu_adam_state(learner.optimizer),
        :config => learner.config,
        :schedule => learner.schedule,
        :global_step => learner.global_step,
        :train_steps => learner.train_steps,
        :last_loss => learner.last_loss,
        :mean_training_loss => mean_training_loss(learner),
        :training_loss_sum => learner.loss_sum,
        :training_loss_count => learner.loss_count,
        :device => training_device_name(learner.device),
        :action_table => copy(action_table),
        :task_metadata => task_metadata,
        :manifest => _checkpoint_manifest_payload(manifest),
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
        :device => training_device_name(learner.device),
        :action_table => copy(PAPER_ACTIONS_MPS),
        :manifest => _checkpoint_manifest_payload(manifest),
    )
    serialize(path, payload)
    return path
end

load_checkpoint(path::AbstractString) = deserialize(path)
