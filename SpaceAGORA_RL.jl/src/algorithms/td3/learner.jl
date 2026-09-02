Base.@kwdef struct TD3Config
    actor_learning_rate::Float64 = 3e-4
    critic_learning_rate::Float64 = 3e-4
    discount::Float64 = 0.95
    batch_size::Int = 256
    train_frequency::Int = 1
    updates_per_step::Int = 1
    train_start::Int = 25_000
    random_steps::Int = 25_000
    replay_size::Int = 1_000_000
    exploration_noise::Float64 = 0.1
    target_policy_noise::Float64 = 0.2
    target_noise_clip::Float64 = 0.5
    policy_delay::Int = 2
    tau::Float64 = 0.005
    gradient_clip_norm::Float64 = 10.0
    adam_epsilon::Float64 = 1e-6
    adam_beta1::Float64 = 0.9
    adam_beta2::Float64 = 0.999
    hidden_dim::Int = 1024
    obs_dim::Int = 9
    action_dim::Int = 1
    bootstrap_truncated::Bool = true
end

mutable struct TD3Learner
    actor::QNetwork
    target_actor::QNetwork
    critic1::QNetwork
    critic2::QNetwork
    target_critic1::QNetwork
    target_critic2::QNetwork
    replay::ContinuousReplayBuffer
    config::TD3Config
    actor_optimizer::AdamState
    critic1_optimizer::AdamState
    critic2_optimizer::AdamState
    global_step::Int
    train_steps::Int
    actor_updates::Int
    last_actor_loss::Float64
    last_critic1_loss::Float64
    last_critic2_loss::Float64
    last_actor_gradient_norm::Float64
    last_critic1_gradient_norm::Float64
    last_critic2_gradient_norm::Float64
    critic_loss_sum::Float64
    critic_loss_count::Int
    device::AbstractTrainingDevice
end

function _td3_optimizer_config(learning_rate::Real, config::TD3Config)
    return (
        learning_rate=Float64(learning_rate),
        adam_beta1=config.adam_beta1,
        adam_beta2=config.adam_beta2,
        adam_epsilon=config.adam_epsilon,
    )
end

function TD3Learner(rng::AbstractRNG, config::TD3Config=TD3Config();
                    device::AbstractTrainingDevice=CPUTrainingDevice())
    actor = network_to_device(device, init_q_network(
        rng;
        input_dim=config.obs_dim,
        hidden_dim=config.hidden_dim,
        output_dim=config.action_dim,
        output_gain=0.01,
    ))
    critic_input_dim = config.obs_dim + config.action_dim
    critic1 = network_to_device(device, init_q_network(
        rng;
        input_dim=critic_input_dim,
        hidden_dim=config.hidden_dim,
        output_dim=1,
        output_gain=1.0,
    ))
    critic2 = network_to_device(device, init_q_network(
        rng;
        input_dim=critic_input_dim,
        hidden_dim=config.hidden_dim,
        output_dim=1,
        output_gain=1.0,
    ))
    return TD3Learner(
        actor,
        copy(actor),
        critic1,
        critic2,
        copy(critic1),
        copy(critic2),
        ContinuousReplayBuffer(config.obs_dim, config.action_dim, config.replay_size),
        config,
        AdamState(actor, _td3_optimizer_config(config.actor_learning_rate, config)),
        AdamState(critic1, _td3_optimizer_config(config.critic_learning_rate, config)),
        AdamState(critic2, _td3_optimizer_config(config.critic_learning_rate, config)),
        0,
        0,
        0,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        0.0,
        0,
        device,
    )
end

td3_actor_output(actor::QNetwork, observations::AbstractMatrix{<:Real}) =
    tanh.(predict_q(actor, observations))

td3_actor_output(actor::QNetwork, observation::AbstractVector{<:Real}) =
    vec(td3_actor_output(actor, reshape(_as_float32_array(observation), :, 1)))

function td3_normalized_action(learner::TD3Learner, observation::AbstractVector{<:Real})
    observations = to_device_array(
        learner.device,
        reshape(_as_float32_array(observation), :, 1),
    )
    return vec(Float32.(to_cpu_array(td3_actor_output(learner.actor, observations))))
end

function select_action(learner::TD3Learner, observation::AbstractVector{<:Real};
                       rng::AbstractRNG=Random.default_rng(), test::Bool=false)
    normalized = if !test && learner.global_step < learner.config.random_steps
        2f0 .* rand(rng, Float32, learner.config.action_dim) .- 1f0
    else
        td3_normalized_action(learner, observation)
    end
    if !test && learner.global_step >= learner.config.random_steps
        normalized .+= Float32(learner.config.exploration_noise) .*
                      randn(rng, Float32, learner.config.action_dim)
    end
    normalized = clamp.(normalized, -1f0, 1f0)
    learner.config.action_dim == 1 || throw(ArgumentError(
        "the aerobraking TD3 policy requires action_dim=1",
    ))
    return _continuous_action_from_normalized(only(normalized))
end

function observe!(learner::TD3Learner, transition::ContinuousTransition)
    push!(learner.replay, transition)
    return learner
end

function td3_targets(learner::TD3Learner, batch::ContinuousReplayBatch,
                     rng::AbstractRNG; target_noise=nothing)
    next_observations = to_device_array(learner.device, batch.next_observations)
    target_actions = td3_actor_output(learner.target_actor, next_observations)
    noise_cpu = if target_noise === nothing
        Float32(learner.config.target_policy_noise) .*
        randn(rng, Float32, size(batch.actions))
    else
        size(target_noise) == size(batch.actions) || throw(DimensionMismatch(
            "TD3 target noise must match the replay action batch",
        ))
        Float32.(target_noise)
    end
    noise_cpu = clamp.(noise_cpu,
                       -Float32(learner.config.target_noise_clip),
                       Float32(learner.config.target_noise_clip))
    smoothed_actions = clamp.(
        target_actions .+ to_device_array(learner.device, noise_cpu),
        -1f0,
        1f0,
    )
    critic_inputs = vcat(next_observations, smoothed_actions)
    q1 = vec(to_cpu_array(predict_q(learner.target_critic1, critic_inputs)))
    q2 = vec(to_cpu_array(predict_q(learner.target_critic2, critic_inputs)))
    targets = Vector{Float32}(undef, length(batch.rewards))
    for index in eachindex(targets)
        done = batch.terminated[index] ||
               (!learner.config.bootstrap_truncated && batch.truncated[index])
        bootstrap = done ? 0f0 : min(q1[index], q2[index])
        targets[index] = batch.rewards[index] +
                         Float32(learner.config.discount) * bootstrap
    end
    return targets
end

function _td3_critic_loss_and_gradients(critic::QNetwork,
                                        observations,
                                        normalized_actions,
                                        targets,
                                        device::AbstractTrainingDevice)
    inputs = vcat(observations, normalized_actions)
    predictions = predict_q(critic, inputs)
    target_row = reshape(to_device_array(device, targets), 1, :)
    errors = predictions .- target_row
    loss = Float64(cpu_scalar(mean(abs2, errors)))
    output_delta = (2f0 / Float32(length(targets))) .* errors
    gradients = network_gradients_from_output_delta(critic, inputs, output_delta)
    return loss, gradients
end

function _td3_actor_loss_and_gradients(learner::TD3Learner, observations)
    raw_actions = predict_q(learner.actor, observations)
    normalized_actions = tanh.(raw_actions)
    critic_inputs = vcat(observations, normalized_actions)
    q_values = predict_q(learner.critic1, critic_inputs)
    batch_size = size(observations, 2)
    critic_output_delta = fill!(similar(q_values), -1f0 / Float32(batch_size))
    _, critic_input_delta = network_gradients_and_input_delta(
        learner.critic1,
        critic_inputs,
        critic_output_delta,
    )
    action_delta = @view critic_input_delta[(learner.config.obs_dim + 1):end, :]
    actor_output_delta = action_delta .* (1f0 .- normalized_actions .* normalized_actions)
    actor_gradients = network_gradients_from_output_delta(
        learner.actor,
        observations,
        actor_output_delta,
    )
    return -Float64(cpu_scalar(mean(q_values))), actor_gradients
end

function _td3_clip_gradients!(gradients::QNetworkGradients, limit::Real)
    norm = gradient_norm(gradients)
    if isfinite(norm) && norm > limit
        scale_gradients!(gradients, limit / norm)
    end
    return norm
end

function train_step!(learner::TD3Learner, rng::AbstractRNG)
    batch = sample_batch(learner.replay, learner.config.batch_size, rng)
    observations = to_device_array(learner.device, batch.observations)
    action_scale = Float32((CONTINUOUS_ACTION_HIGH_MPS - CONTINUOUS_ACTION_LOW_MPS) / 2)
    action_midpoint = Float32((CONTINUOUS_ACTION_HIGH_MPS + CONTINUOUS_ACTION_LOW_MPS) / 2)
    normalized_actions_cpu = clamp.((batch.actions .- action_midpoint) ./ action_scale,
                                    -1f0, 1f0)
    normalized_actions = to_device_array(learner.device, normalized_actions_cpu)
    targets = td3_targets(learner, batch, rng)
    critic1_loss, critic1_gradients = _td3_critic_loss_and_gradients(
        learner.critic1,
        observations,
        normalized_actions,
        targets,
        learner.device,
    )
    critic2_loss, critic2_gradients = _td3_critic_loss_and_gradients(
        learner.critic2,
        observations,
        normalized_actions,
        targets,
        learner.device,
    )
    learner.last_critic1_gradient_norm = _td3_clip_gradients!(
        critic1_gradients,
        learner.config.gradient_clip_norm,
    )
    learner.last_critic2_gradient_norm = _td3_clip_gradients!(
        critic2_gradients,
        learner.config.gradient_clip_norm,
    )
    adam_update!(learner.critic1, critic1_gradients, learner.critic1_optimizer)
    adam_update!(learner.critic2, critic2_gradients, learner.critic2_optimizer)
    learner.train_steps += 1
    learner.last_critic1_loss = critic1_loss
    learner.last_critic2_loss = critic2_loss
    learner.critic_loss_sum += (critic1_loss + critic2_loss) / 2
    learner.critic_loss_count += 1

    if learner.train_steps % learner.config.policy_delay == 0
        actor_loss, actor_gradients = _td3_actor_loss_and_gradients(learner, observations)
        learner.last_actor_gradient_norm = _td3_clip_gradients!(
            actor_gradients,
            learner.config.gradient_clip_norm,
        )
        adam_update!(learner.actor, actor_gradients, learner.actor_optimizer)
        polyak_update!(learner.target_actor, learner.actor, learner.config.tau)
        polyak_update!(learner.target_critic1, learner.critic1, learner.config.tau)
        polyak_update!(learner.target_critic2, learner.critic2, learner.config.tau)
        learner.actor_updates += 1
        learner.last_actor_loss = actor_loss
    end
    return (critic1_loss=critic1_loss, critic2_loss=critic2_loss,
            actor_loss=learner.last_actor_loss)
end

mean_training_loss(learner::TD3Learner) = learner.critic_loss_count == 0 ? NaN :
    learner.critic_loss_sum / learner.critic_loss_count

function maybe_train!(learner::TD3Learner, rng::AbstractRNG)
    if length(learner.replay) >= learner.config.batch_size &&
       learner.global_step >= learner.config.train_start &&
       learner.global_step % learner.config.train_frequency == 0
        result = nothing
        for _ in 1:learner.config.updates_per_step
            result = train_step!(learner, rng)
        end
        return result
    end
    return nothing
end
