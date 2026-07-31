Base.@kwdef struct A2CConfig
    learning_rate::Float64 = 1e-4
    discount::Float64 = 0.95
    segment_length::Int = 10
    train_start::Int = 0
    entropy_coef::Float64 = 0.1
    value_coef::Float64 = 0.5
    gradient_clip_norm::Float64 = 0.5
    adam_epsilon::Float64 = 1e-6
    adam_beta1::Float64 = 0.9
    adam_beta2::Float64 = 0.99
    hidden_dim::Int = 1024
    obs_dim::Int = 9
    action_dim::Int = 13
end

mutable struct A2CLearner
    actor::QNetwork
    critic::QNetwork
    config::A2CConfig
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
    device::AbstractTrainingDevice
end

function A2CLearner(rng::AbstractRNG, config::A2CConfig=A2CConfig();
                    device::AbstractTrainingDevice=CPUTrainingDevice())
    actor = network_to_device(device, init_q_network(rng; input_dim=config.obs_dim,
                                                     hidden_dim=config.hidden_dim,
                                                     output_dim=config.action_dim))
    critic = network_to_device(device, init_q_network(rng; input_dim=config.obs_dim,
                                                      hidden_dim=config.hidden_dim,
                                                      output_dim=1))
    return A2CLearner(
        actor,
        critic,
        config,
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
        device,
    )
end

function softmax_columns(logits::AbstractMatrix{Float32})
    shifted = logits .- maximum(logits; dims=1)
    exp_logits = exp.(shifted)
    return exp_logits ./ sum(exp_logits; dims=1)
end

function sample_categorical(probs::AbstractVector{<:Real}, rng::AbstractRNG)
    threshold = rand(rng)
    cumulative = 0.0
    for (idx, prob) in pairs(probs)
        cumulative += Float64(prob)
        threshold <= cumulative && return idx
    end
    return length(probs)
end

function actor_action(actor::QNetwork, observation::AbstractVector{<:Real}, rng::AbstractRNG;
                      test::Bool=false)
    logits = predict_q(actor, observation)
    test && return argmax(logits)
    probs = vec(softmax_columns(reshape(Float32.(logits), :, 1)))
    return sample_categorical(probs, rng)
end

function select_action(learner::A2CLearner, observation::AbstractVector{<:Real};
                       rng::AbstractRNG=Random.default_rng(), test::Bool=false)
    learner.global_step += 1
    return actor_action(cpu_network(learner.actor), observation, rng; test=test)
end

function value_predictions(learner::A2CLearner, observations::AbstractMatrix{Float32})
    obs = to_device_array(learner.device, observations)
    return Float32.(vec(to_cpu_array(predict_q(learner.critic, obs))))
end

function a2c_loss_and_gradients(learner::A2CLearner, batch::A2CRolloutBatch)
    n = length(batch)
    n > 0 || throw(ArgumentError("A2C rollout batch must be nonempty"))
    observations = to_device_array(learner.device, batch.observations)
    returns = to_device_array(learner.device, batch.returns)

    logits = predict_q(learner.actor, observations)
    values = vec(predict_q(learner.critic, observations))
    probs = softmax_columns(logits)
    log_probs = log.(max.(probs, eps(Float32)))
    encoded_actions = to_device_array(learner.device, onehot_actions(batch.actions, learner.config.action_dim))

    raw_advantages = returns .- values
    raw_advantages_cpu = Float32.(to_cpu_array(raw_advantages))
    advantage_scale = 1f0 + Float32(std(raw_advantages_cpu; corrected=false))
    advantages_cpu = raw_advantages_cpu ./ advantage_scale
    advantages = to_device_array(learner.device, advantages_cpu)

    selected_log_probs = vec(sum(log_probs .* encoded_actions; dims=1))
    entropy_cols = vec(-sum(probs .* log_probs; dims=1))
    policy_loss = -mean(to_cpu_array(selected_log_probs) .* advantages_cpu)
    entropy_loss = mean(to_cpu_array(entropy_cols))
    value_loss = mean(abs2, advantages_cpu)
    total_loss = policy_loss - learner.config.entropy_coef * entropy_loss +
                 learner.config.value_coef * value_loss

    n32 = Float32(n)
    advantage_row = reshape(advantages, 1, :)
    entropy_row = reshape(-sum(probs .* log_probs; dims=1), 1, :)
    dlogits = ((probs .- encoded_actions) .* advantage_row) ./ n32
    dlogits .+= (Float32(learner.config.entropy_coef) / n32) .* probs .* (log_probs .+ entropy_row)

    value_delta = (-2f0 * Float32(learner.config.value_coef) / (n32 * advantage_scale)) .*
                  reshape(advantages, 1, :)

    actor_grads = network_gradients_from_output_delta(learner.actor, observations, dlogits)
    critic_grads = network_gradients_from_output_delta(learner.critic, observations, value_delta)
    return total_loss, policy_loss, entropy_loss, value_loss, actor_grads, critic_grads
end

function train_step!(learner::A2CLearner, batch::A2CRolloutBatch)
    loss, policy_loss, entropy_loss, value_loss, actor_grads, critic_grads =
        a2c_loss_and_gradients(learner, batch)

    actor_norm = gradient_norm(actor_grads)
    if isfinite(actor_norm) && actor_norm > learner.config.gradient_clip_norm
        scale_gradients!(actor_grads, learner.config.gradient_clip_norm / actor_norm)
    end
    critic_norm = gradient_norm(critic_grads)
    if isfinite(critic_norm) && critic_norm > learner.config.gradient_clip_norm
        scale_gradients!(critic_grads, learner.config.gradient_clip_norm / critic_norm)
    end

    adam_update!(learner.actor, actor_grads, learner.actor_optimizer)
    adam_update!(learner.critic, critic_grads, learner.critic_optimizer)
    learner.train_steps += 1
    learner.last_loss = loss
    learner.loss_sum += loss
    learner.loss_count += 1
    learner.last_policy_loss = policy_loss
    learner.last_entropy = entropy_loss
    learner.last_value_loss = value_loss
    return loss
end

mean_training_loss(learner::A2CLearner) =
    learner.loss_count == 0 ? NaN : learner.loss_sum / learner.loss_count

function maybe_train!(learner::A2CLearner, batch::A2CRolloutBatch)
    learner.global_step >= learner.config.train_start || return nothing
    return train_step!(learner, batch)
end
