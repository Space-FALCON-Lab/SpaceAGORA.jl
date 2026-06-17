Base.@kwdef struct AADSHeuristicPolicy <: AbstractPolicy
    target_heat_rate_w_cm2::Float64 = 0.15
end

function predicted_heat_rate(config, state, action_index::Integer)
    action = action_from_index(action_index)
    next_hp = clamp(apply_apoapsis_maneuver(config, state, action), 50e3, 180e3)
    rp = config.mars_radius_m + next_hp
    a = (state.apoapsis_radius_m + rp) / 2
    v = sqrt(config.mu_m3_s2 * (2 / rp - 1 / a))
    return paper_heat_rate_w_cm2(config, next_hp, v)
end

function policy_action_index(policy::AADSHeuristicPolicy, config, state, observation::PaperObservation, rng::AbstractRNG)
    zero_idx = zero_action_index()
    zero_heat = predicted_heat_rate(config, state, zero_idx)
    if config.reward_config.heat_low_w_cm2 <= zero_heat <= config.reward_config.heat_high_w_cm2
        return zero_idx
    end

    candidates = zero_heat < config.reward_config.heat_low_w_cm2 ?
                 findall(<(0.0), PAPER_ACTIONS_MPS) :
                 findall(>(0.0), PAPER_ACTIONS_MPS)
    isempty(candidates) && return zero_idx

    _, local_idx = findmin(abs.(predicted_heat_rate.(Ref(config), Ref(state), candidates) .-
                                policy.target_heat_rate_w_cm2))
    return candidates[local_idx]
end
