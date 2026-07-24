Base.@kwdef struct DDQNConfig
    learning_rate::Float64 = 1e-4
    discount::Float64 = 0.95
    batch_size::Int = 256
    train_frequency::Int = 5
    train_start::Int = 10_000
    replay_size::Int = 1_000_000
    target_update::Int = 10_000
    gradient_clip_norm::Float64 = 10.0
    adam_epsilon::Float64 = 1e-6
    adam_beta1::Float64 = 0.9
    adam_beta2::Float64 = 0.999
    hidden_dim::Int = 1024
    obs_dim::Int = 9
    action_dim::Int = 13
    bootstrap_truncated::Bool = true
end

mutable struct AdamState
    m::QNetwork
    v::QNetwork
    t::Int
    beta1::Float32
    beta2::Float32
    epsilon::Float32
    learning_rate::Float32
end

function AdamState(net::QNetwork, config)
    return AdamState(zero_network_like(net), zero_network_like(net), 0,
                     Float32(config.adam_beta1), Float32(config.adam_beta2),
                     Float32(config.adam_epsilon), Float32(config.learning_rate))
end

mutable struct DDQNLearner
    online::QNetwork
    target::QNetwork
    replay::ReplayBuffer
    schedule::EpsilonSchedule
    config::DDQNConfig
    optimizer::AdamState
    global_step::Int
    train_steps::Int
    last_loss::Float64
    device::AbstractTrainingDevice
end

function DDQNLearner(rng::AbstractRNG, config::DDQNConfig=DDQNConfig();
                     schedule::EpsilonSchedule=EpsilonSchedule(decay_start_step=config.train_start),
                     device::AbstractTrainingDevice=CPUTrainingDevice())
    online = network_to_device(device, init_q_network(rng; input_dim=config.obs_dim,
                                                      hidden_dim=config.hidden_dim,
                                                      output_dim=config.action_dim))
    target = copy(online)
    replay = ReplayBuffer(config.obs_dim, config.replay_size)
    optimizer = AdamState(online, config)
    return DDQNLearner(online, target, replay, schedule, config, optimizer, 0, 0, NaN, device)
end

function selected_q_targets(online_q_next::AbstractMatrix{<:Real},
                            target_q_next::AbstractMatrix{<:Real},
                            rewards::AbstractVector{<:Real},
                            terminated::AbstractVector{Bool},
                            truncated::AbstractVector{Bool},
                            gamma::Real;
                            bootstrap_truncated::Bool=true)
    batch_size = length(rewards)
    targets = Vector{Float32}(undef, batch_size)
    for i in 1:batch_size
        online_action = argmax(view(online_q_next, :, i))
        done = terminated[i] || (!bootstrap_truncated && truncated[i])
        bootstrap = done ? 0.0 : Float64(target_q_next[online_action, i])
        targets[i] = Float32(Float64(rewards[i]) + Float64(gamma) * bootstrap)
    end
    return targets
end

compute_ddqn_targets(args...; kwargs...) = selected_q_targets(args...; kwargs...)

function greedy_action_index(learner::DDQNLearner, observation::AbstractVector{<:Real})
    obs = reshape(_as_float32_array(observation), :, 1)
    q = to_cpu_array(predict_q(learner.online, to_device_array(learner.device, obs)))
    return argmax(view(q, :, 1))
end

function select_action(learner::DDQNLearner, observation::AbstractVector{<:Real};
                       rng::AbstractRNG=Random.default_rng(), test::Bool=false)
    learner.global_step += 1
    eps = test ? 0.0 : epsilon_value(learner.schedule, learner.global_step)
    if !test && rand(rng) < eps
        return rand(rng, 1:learner.config.action_dim)
    end
    return greedy_action_index(learner, observation)
end

function observe!(learner::DDQNLearner, transition::Transition)
    push!(learner.replay, transition)
    return learner
end

function adam_update_array!(param, grad, m, v, opt::AdamState)
    opt.t >= 1 || throw(ArgumentError("Adam step counter must be incremented before update"))
    b1 = opt.beta1
    b2 = opt.beta2
    m .= b1 .* m .+ (1f0 - b1) .* grad
    v .= b2 .* v .+ (1f0 - b2) .* (grad .* grad)
    one_minus_b1t = 1f0 - b1^opt.t
    one_minus_b2t = 1f0 - b2^opt.t
    param .-= opt.learning_rate .* (m ./ one_minus_b1t) ./
              (sqrt.(v ./ one_minus_b2t) .+ opt.epsilon)
    return param
end

function adam_update!(net::QNetwork, grads::QNetworkGradients, opt::AdamState)
    opt.t += 1
    adam_update_array!(net.W1, grads.W1, opt.m.W1, opt.v.W1, opt)
    adam_update_array!(net.b1, grads.b1, opt.m.b1, opt.v.b1, opt)
    adam_update_array!(net.W2, grads.W2, opt.m.W2, opt.v.W2, opt)
    adam_update_array!(net.b2, grads.b2, opt.m.b2, opt.v.b2, opt)
    adam_update_array!(net.W3, grads.W3, opt.m.W3, opt.v.W3, opt)
    adam_update_array!(net.b3, grads.b3, opt.m.b3, opt.v.b3, opt)
    return net
end

function train_step!(learner::DDQNLearner, rng::AbstractRNG)
    batch = sample_batch(learner.replay, learner.config.batch_size, rng)
    observations = to_device_array(learner.device, batch.observations)
    next_observations = to_device_array(learner.device, batch.next_observations)
    online_next = predict_q(learner.online, next_observations)
    target_next = predict_q(learner.target, next_observations)
    targets = compute_ddqn_targets(to_cpu_array(online_next), to_cpu_array(target_next),
                                   batch.rewards, batch.terminated, batch.truncated,
                                   learner.config.discount;
                                   bootstrap_truncated=learner.config.bootstrap_truncated)
    loss, grads = network_loss_and_gradients(learner.online, observations, batch.actions, targets;
                                             device=learner.device)
    norm = gradient_norm(grads)
    if isfinite(norm) && norm > learner.config.gradient_clip_norm
        scale_gradients!(grads, learner.config.gradient_clip_norm / norm)
    end
    adam_update!(learner.online, grads, learner.optimizer)
    learner.train_steps += 1
    learner.last_loss = loss
    maybe_update_target!(learner)
    return loss
end

function maybe_train!(learner::DDQNLearner, rng::AbstractRNG)
    if length(learner.replay) >= learner.config.batch_size &&
       learner.global_step >= learner.config.train_start &&
       learner.global_step % learner.config.train_frequency == 0
        return train_step!(learner, rng)
    end
    return nothing
end
