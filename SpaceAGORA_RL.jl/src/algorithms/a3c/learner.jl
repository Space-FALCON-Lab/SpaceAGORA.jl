Base.@kwdef struct A3CConfig
    learning_rate::Float64 = 1e-4
    discount::Float64 = 0.95
    t_max::Int = 10
    entropy_coef::Float64 = 0.01
    value_coef::Float64 = 0.5
    normalize_advantages::Bool = false
    gradient_clip_norm::Float64 = 0.5
    adam_epsilon::Float64 = 1e-6
    adam_beta1::Float64 = 0.9
    adam_beta2::Float64 = 0.99
    hidden_dim::Int = 1024
    obs_dim::Int = 9
    action_dim::Int = 13
    max_policy_lag::Int = -1
end

mutable struct A3CLearner
    actor::QNetwork
    critic::QNetwork
    config::A3CConfig
    action_config::ActorCriticActionConfig
    actor_optimizer::AdamState
    critic_optimizer::AdamState
    global_step::Int
    train_steps::Int
    last_loss::Float64
    loss_sum::Float64
    loss_count::Int
    last_policy_loss::Float64
    last_entropy::Float64
    last_value_loss::Float64
    last_explained_variance::Float64
    last_actor_gradient_norm::Float64
    last_critic_gradient_norm::Float64
    policy_version::Int
    gradient_lag_sum::Int
    gradient_lag_count::Int
    max_observed_policy_lag::Int
    dropped_stale_updates::Int
    device::AbstractTrainingDevice
end

mutable struct A3CLocalModel
    actor::QNetwork
    critic::QNetwork
    policy_version::Int
    action_config::ActorCriticActionConfig
end

function A3CLearner(rng::AbstractRNG, config::A3CConfig=A3CConfig();
                    device::AbstractTrainingDevice=CPUTrainingDevice(),
                    action_config::ActorCriticActionConfig=ActorCriticActionConfig())
    actor = network_to_device(device, init_actor_network(
        rng,
        config.obs_dim,
        config.hidden_dim,
        config.action_dim,
        action_config,
    ))
    critic = network_to_device(device, init_q_network(
        rng;
        input_dim=config.obs_dim,
        hidden_dim=config.hidden_dim,
        output_dim=1,
        output_gain=1.0,
    ))
    return A3CLearner(
        actor,
        critic,
        config,
        action_config,
        AdamState(actor, config),
        AdamState(critic, config),
        0,
        0,
        NaN,
        0.0,
        0,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        0,
        0,
        0,
        0,
        0,
        device,
    )
end

A3CLocalModel(learner::A3CLearner) =
    A3CLocalModel(
        cpu_network(learner.actor),
        cpu_network(learner.critic),
        learner.policy_version,
        learner.action_config,
    )

function sync_local_model!(local_model::A3CLocalModel, learner::A3CLearner)
    local_model.actor = cpu_network(learner.actor)
    local_model.critic = cpu_network(learner.critic)
    local_model.policy_version = learner.policy_version
    return local_model
end

function select_action(local_model::A3CLocalModel, observation::AbstractVector{<:Real};
                       rng::AbstractRNG=Random.default_rng(), test::Bool=false)
    return actor_critic_action(
        local_model.actor,
        local_model.action_config,
        observation,
        rng;
        test=test,
    )
end

function a3c_rollout_batch(rollout::A3CRollout, local_model::A3CLocalModel,
                           config::A3CConfig)
    isempty(rollout) && throw(ArgumentError("A3C rollout must be nonempty"))
    transitions = rollout.transitions
    n = length(transitions)
    observations = Matrix{Float32}(undef, config.obs_dim, n)
    actions = Vector{Int}(undef, n)
    returns = Vector{Float32}(undef, n)
    last_transition = last(transitions)
    running = last_transition.terminated ?
              0f0 :
              Float32(only(predict_q(local_model.critic, last_transition.next_observation)))
    gamma = Float32(config.discount)
    for index in n:-1:1
        transition = transitions[index]
        if index < n && transition.terminated
            throw(ArgumentError("an A3C rollout cannot continue after a terminal transition"))
        end
        observations[:, index] .= transition.observation
        actions[index] = transition.action_index
        running = transition.reward + gamma * running
        returns[index] = running
    end
    if uses_continuous_actions(local_model.action_config)
        length(rollout.actions_mps) == n || throw(DimensionMismatch(
            "A3C continuous-action rollout values must match its transitions",
        ))
        return ContinuousA3CRolloutBatch(
            observations,
            copy(rollout.actions_mps),
            returns,
            local_model.policy_version,
        )
    end
    return A3CRolloutBatch(observations, actions, returns, local_model.policy_version)
end

function a3c_loss_and_gradients(local_model::A3CLocalModel, config::A3CConfig,
                                batch::A3CRolloutBatch)
    n = length(batch)
    n > 0 || throw(ArgumentError("A3C rollout batch must be nonempty"))
    batch.policy_version == local_model.policy_version || throw(ArgumentError(
        "A3C rollout policy version $(batch.policy_version) does not match local version $(local_model.policy_version)",
    ))
    logits = predict_q(local_model.actor, batch.observations)
    values = vec(predict_q(local_model.critic, batch.observations))
    probs = softmax_columns(logits)
    log_probs = log.(max.(probs, eps(Float32)))
    encoded_actions = onehot_actions(batch.actions, config.action_dim)
    raw_advantages = batch.returns .- values
    advantages = if config.normalize_advantages && n > 1
        advantage_mean = Float32(mean(raw_advantages))
        advantage_scale = Float32(std(raw_advantages; corrected=false)) + eps(Float32)
        (raw_advantages .- advantage_mean) ./ advantage_scale
    else
        copy(raw_advantages)
    end

    selected_log_probs = vec(sum(log_probs .* encoded_actions; dims=1))
    entropy_cols = vec(-sum(probs .* log_probs; dims=1))
    policy_loss = Float64(-mean(selected_log_probs .* advantages))
    entropy = Float64(mean(entropy_cols))
    value_loss = Float64(mean(abs2, raw_advantages))
    total_loss = policy_loss - config.entropy_coef * entropy + config.value_coef * value_loss

    n32 = Float32(n)
    advantage_row = reshape(advantages, 1, :)
    entropy_row = reshape(-sum(probs .* log_probs; dims=1), 1, :)
    dlogits = ((probs .- encoded_actions) .* advantage_row) ./ n32
    dlogits .+= (Float32(config.entropy_coef) / n32) .* probs .* (log_probs .+ entropy_row)
    value_delta = (-2f0 * Float32(config.value_coef) / n32) .* reshape(raw_advantages, 1, :)

    actor_grads = network_gradients_from_output_delta(local_model.actor, batch.observations, dlogits)
    critic_grads = network_gradients_from_output_delta(local_model.critic, batch.observations, value_delta)
    return (
        total_loss,
        policy_loss,
        entropy,
        value_loss,
        Float64(var(batch.returns; corrected=false)),
        Float64(var(raw_advantages; corrected=false)),
        actor_grads,
        critic_grads,
    )
end

function a3c_loss_and_gradients(local_model::A3CLocalModel, config::A3CConfig,
                                batch::ContinuousA3CRolloutBatch)
    uses_continuous_actions(local_model.action_config) || throw(ArgumentError(
        "continuous A3C batches require a continuous-action local model",
    ))
    n = length(batch)
    n > 0 || throw(ArgumentError("A3C rollout batch must be nonempty"))
    batch.policy_version == local_model.policy_version || throw(ArgumentError(
        "A3C rollout policy version $(batch.policy_version) does not match local version $(local_model.policy_version)",
    ))
    values = vec(predict_q(local_model.critic, batch.observations))
    raw_advantages = batch.returns .- values
    advantages = if config.normalize_advantages && n > 1
        advantage_mean = Float32(mean(raw_advantages))
        advantage_scale = Float32(std(raw_advantages; corrected=false)) + eps(Float32)
        (raw_advantages .- advantage_mean) ./ advantage_scale
    else
        copy(raw_advantages)
    end
    policy_loss, entropy, actor_delta = continuous_policy_loss_and_output_delta(
        local_model.actor,
        local_model.action_config,
        batch.observations,
        batch.actions_mps,
        advantages,
        config.entropy_coef,
    )
    value_loss = Float64(mean(abs2, raw_advantages))
    total_loss = policy_loss - config.entropy_coef * entropy + config.value_coef * value_loss
    value_delta = (-2f0 * Float32(config.value_coef) / Float32(n)) .*
                  reshape(raw_advantages, 1, :)
    actor_grads = network_gradients_from_output_delta(
        local_model.actor,
        batch.observations,
        actor_delta,
    )
    critic_grads = network_gradients_from_output_delta(
        local_model.critic,
        batch.observations,
        value_delta,
    )
    return (
        total_loss,
        policy_loss,
        entropy,
        value_loss,
        Float64(var(batch.returns; corrected=false)),
        Float64(var(raw_advantages; corrected=false)),
        actor_grads,
        critic_grads,
    )
end

function _a3c_gradients_to_device(device::AbstractTrainingDevice,
                                  gradients::QNetworkGradients)
    return QNetworkGradients(
        to_device_array(device, gradients.W1),
        to_device_array(device, gradients.b1),
        to_device_array(device, gradients.W2),
        to_device_array(device, gradients.b2),
        to_device_array(device, gradients.W3),
        to_device_array(device, gradients.b3),
    )
end

function a3c_update!(learner::A3CLearner, local_model::A3CLocalModel,
                     batch::Union{A3CRolloutBatch,ContinuousA3CRolloutBatch})
    loss, policy_loss, entropy, value_loss, return_variance, residual_variance,
        actor_grads, critic_grads = a3c_loss_and_gradients(local_model, learner.config, batch)
    lag = learner.policy_version - local_model.policy_version
    lag >= 0 || throw(ArgumentError("A3C local policy version cannot exceed the global version"))
    learner.gradient_lag_sum += lag
    learner.gradient_lag_count += 1
    learner.max_observed_policy_lag = max(learner.max_observed_policy_lag, lag)

    if learner.config.max_policy_lag >= 0 && lag > learner.config.max_policy_lag
        learner.dropped_stale_updates += 1
        sync_local_model!(local_model, learner)
        return nothing
    end

    actor_norm = gradient_norm(actor_grads)
    if isfinite(actor_norm) && actor_norm > learner.config.gradient_clip_norm
        scale_gradients!(actor_grads, learner.config.gradient_clip_norm / actor_norm)
    end
    critic_norm = gradient_norm(critic_grads)
    if isfinite(critic_norm) && critic_norm > learner.config.gradient_clip_norm
        scale_gradients!(critic_grads, learner.config.gradient_clip_norm / critic_norm)
    end

    adam_update!(
        learner.actor,
        _a3c_gradients_to_device(learner.device, actor_grads),
        learner.actor_optimizer,
    )
    adam_update!(
        learner.critic,
        _a3c_gradients_to_device(learner.device, critic_grads),
        learner.critic_optimizer,
    )
    learner.train_steps += 1
    learner.policy_version += 1
    learner.last_loss = loss
    learner.loss_sum += loss
    learner.loss_count += 1
    learner.last_policy_loss = policy_loss
    learner.last_entropy = entropy
    learner.last_value_loss = value_loss
    learner.last_explained_variance = return_variance > eps(Float32) ?
                                      1.0 - residual_variance / return_variance : NaN
    learner.last_actor_gradient_norm = actor_norm
    learner.last_critic_gradient_norm = critic_norm
    sync_local_model!(local_model, learner)
    return loss
end

mean_training_loss(learner::A3CLearner) =
    learner.loss_count == 0 ? NaN : learner.loss_sum / learner.loss_count

mean_policy_lag(learner::A3CLearner) = learner.gradient_lag_count == 0 ?
                                       NaN :
                                       learner.gradient_lag_sum / learner.gradient_lag_count
