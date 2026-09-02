const ACTOR_CRITIC_ACTION_MODES = (:discrete, :continuous)
const CONTINUOUS_ACTION_LOW_MPS = minimum(PAPER_ACTIONS_MPS)
const CONTINUOUS_ACTION_HIGH_MPS = maximum(PAPER_ACTIONS_MPS)

Base.@kwdef struct ActorCriticActionConfig
    mode::Symbol = :discrete
    initial_log_std::Float64 = -1.0
    log_std_min::Float64 = -5.0
    log_std_max::Float64 = 1.0
end

canonical_actor_critic_action_mode(value) =
    Symbol(replace(lowercase(strip(String(value))), "-" => "_"))

uses_continuous_actions(config::ActorCriticActionConfig) = config.mode == :continuous

function _continuous_action_from_normalized(value::Real)
    normalized = clamp(Float64(value), -1.0, 1.0)
    midpoint = (CONTINUOUS_ACTION_LOW_MPS + CONTINUOUS_ACTION_HIGH_MPS) / 2
    half_range = (CONTINUOUS_ACTION_HIGH_MPS - CONTINUOUS_ACTION_LOW_MPS) / 2
    return continuous_action_from_delta_v(midpoint + half_range * normalized)
end

function _normalized_continuous_action(delta_v_mps::Real)
    midpoint = (CONTINUOUS_ACTION_LOW_MPS + CONTINUOUS_ACTION_HIGH_MPS) / 2
    half_range = (CONTINUOUS_ACTION_HIGH_MPS - CONTINUOUS_ACTION_LOW_MPS) / 2
    return Float32(clamp((Float64(delta_v_mps) - midpoint) / half_range, -1.0, 1.0))
end

function init_actor_network(rng::AbstractRNG, input_dim::Int, hidden_dim::Int,
                            action_dim::Int, action_config::ActorCriticActionConfig)
    output_dim = uses_continuous_actions(action_config) ? 2 : action_dim
    actor = init_q_network(
        rng;
        input_dim=input_dim,
        hidden_dim=hidden_dim,
        output_dim=output_dim,
        output_gain=0.01,
    )
    if uses_continuous_actions(action_config)
        actor.b3[2] = Float32(action_config.initial_log_std)
    end
    return actor
end

function continuous_actor_action(actor::QNetwork, action_config::ActorCriticActionConfig,
                                 observation::AbstractVector{<:Real}, rng::AbstractRNG;
                                 test::Bool=false)
    parameters = predict_q(actor, observation)
    length(parameters) == 2 || throw(DimensionMismatch(
        "continuous actor must output mean and log standard deviation",
    ))
    mean_value = Float64(parameters[1])
    log_std = clamp(
        Float64(parameters[2]),
        action_config.log_std_min,
        action_config.log_std_max,
    )
    latent = test ? mean_value : mean_value + exp(log_std) * randn(rng)
    return _continuous_action_from_normalized(tanh(latent))
end

function actor_critic_action(actor::QNetwork, action_config::ActorCriticActionConfig,
                             observation::AbstractVector{<:Real}, rng::AbstractRNG;
                             test::Bool=false)
    if uses_continuous_actions(action_config)
        return continuous_actor_action(actor, action_config, observation, rng; test=test)
    end
    return actor_action(actor, observation, rng; test=test)
end

function _actor_critic_action_command(selected;
                                      protected::Bool=false,
                                      policy_version::Int=0)
    if selected isa AerobrakingAction
        return (action=selected, protected=protected, policy_version=policy_version)
    end
    return (action_index=Int(selected), protected=protected, policy_version=policy_version)
end

function continuous_policy_loss_and_output_delta(
    actor::QNetwork,
    action_config::ActorCriticActionConfig,
    observations::AbstractMatrix{Float32},
    actions_mps::AbstractVector{Float32},
    advantages::AbstractVector{Float32},
    entropy_coef::Real,
)
    n = size(observations, 2)
    length(actions_mps) == n ||
        throw(DimensionMismatch("continuous actions length must match batch size"))
    length(advantages) == n ||
        throw(DimensionMismatch("advantages length must match batch size"))
    parameters = predict_q(actor, observations)
    size(parameters) == (2, n) || throw(DimensionMismatch(
        "continuous actor must output a 2 by batch-size parameter matrix",
    ))

    raw_log_std = vec(@view parameters[2, :])
    log_std = clamp.(raw_log_std,
                     Float32(action_config.log_std_min),
                     Float32(action_config.log_std_max))
    means = vec(@view parameters[1, :])
    bounded_actions = clamp.(
        actions_mps,
        Float32(CONTINUOUS_ACTION_LOW_MPS),
        Float32(CONTINUOUS_ACTION_HIGH_MPS),
    )
    action_midpoint = Float32((CONTINUOUS_ACTION_LOW_MPS + CONTINUOUS_ACTION_HIGH_MPS) / 2)
    action_half_range = Float32((CONTINUOUS_ACTION_HIGH_MPS - CONTINUOUS_ACTION_LOW_MPS) / 2)
    normalized_actions = (bounded_actions .- action_midpoint) ./ action_half_range
    stable_actions = clamp.(normalized_actions, -1f0 + 1f-6, 1f0 - 1f-6)
    latent_actions = atanh.(stable_actions)
    centered = latent_actions .- means
    inverse_variance = exp.(-2f0 .* log_std)
    squared_standardized = centered .* centered .* inverse_variance
    action_scale = Float32((CONTINUOUS_ACTION_HIGH_MPS - CONTINUOUS_ACTION_LOW_MPS) / 2)
    log_probabilities = -0.5f0 .* squared_standardized .- log_std .-
                        Float32(0.5 * log(2pi)) .- log(action_scale) .-
                        log.(1f0 .- stable_actions .* stable_actions .+ 1f-6)
    gaussian_entropies = log_std .+ Float32(0.5 * log(2pi * exp(1)))

    policy_loss = Float64(-mean(to_cpu_array(log_probabilities) .* to_cpu_array(advantages)))
    entropy = Float64(mean(to_cpu_array(gaussian_entropies)))
    n32 = Float32(n)
    dmeans = advantages .* (means .- latent_actions) .* inverse_variance ./ n32
    dlog_std = (advantages .* (1f0 .- squared_standardized) .-
                Float32(entropy_coef)) ./ n32
    unclamped = Float32.((raw_log_std .>= Float32(action_config.log_std_min)) .&
                         (raw_log_std .<= Float32(action_config.log_std_max)))
    dlog_std .*= unclamped
    output_delta = vcat(reshape(dmeans, 1, :), reshape(dlog_std, 1, :))
    return policy_loss, entropy, output_delta
end
