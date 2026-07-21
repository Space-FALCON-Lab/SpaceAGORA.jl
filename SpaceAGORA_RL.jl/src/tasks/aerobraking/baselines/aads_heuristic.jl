Base.@kwdef struct AADSHeuristicPolicy <: AbstractPolicy
    target_heat_rate_w_cm2::Float64 = 0.15
    bisection_iterations::Int = 40
end

function predicted_heat_rate(config, state, action_index::Integer)
    return deterministic_predicted_heat_rate(config, state, action_from_index(action_index))
end

function predicted_heat_rate(config, state, action::AerobrakingAction)
    return deterministic_predicted_heat_rate(config, state, action)
end

function _aads_bisection_delta_v(policy::AADSHeuristicPolicy, config, state,
                                 lo::Real, hi::Real)
    lo_f = Float64(lo)
    hi_f = Float64(hi)
    target = policy.target_heat_rate_w_cm2
    f_lo = predicted_heat_rate(config, state, action_from_delta_v(lo_f)) - target
    f_hi = predicted_heat_rate(config, state, action_from_delta_v(hi_f)) - target

    if sign(f_lo) == sign(f_hi)
        return abs(f_lo) <= abs(f_hi) ? lo_f : hi_f
    end

    left = lo_f
    right = hi_f
    f_left = f_lo
    for _ in 1:policy.bisection_iterations
        mid = 0.5 * (left + right)
        f_mid = predicted_heat_rate(config, state, action_from_delta_v(mid)) - target
        abs(f_mid) <= 1e-8 && return mid
        if sign(f_mid) == sign(f_left)
            left = mid
            f_left = f_mid
        else
            right = mid
        end
    end
    return 0.5 * (left + right)
end

function policy_action_index(policy::AADSHeuristicPolicy, config, state, observation::PaperObservation, rng::AbstractRNG)
    zero_idx = zero_action_index()
    zero_heat = predicted_heat_rate(config, state, zero_idx)
    if config.reward_config.heat_low_w_cm2 <= zero_heat <= config.reward_config.heat_high_w_cm2
        return zero_idx
    end

    delta_v = zero_heat < config.reward_config.heat_low_w_cm2 ?
              _aads_bisection_delta_v(policy, config, state, minimum(PAPER_ACTIONS_MPS), 0.0) :
              _aads_bisection_delta_v(policy, config, state, 0.0, maximum(PAPER_ACTIONS_MPS))
    return action_from_delta_v(delta_v)
end
