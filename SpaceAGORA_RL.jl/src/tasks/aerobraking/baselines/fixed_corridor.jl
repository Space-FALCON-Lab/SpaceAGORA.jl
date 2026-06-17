Base.@kwdef struct FixedCorridorPolicy <: AbstractPolicy
    low_action_index::Int = nearest_action_index(-0.2)
    high_action_index::Int = nearest_action_index(0.2)
end

function policy_action_index(policy::FixedCorridorPolicy, config, state, observation::PaperObservation, rng::AbstractRNG)
    heat = observation.max_heat_rate_w_cm2
    target_error = observation.apoapsis_radius_m - config.final_apoapsis_radius_m
    status = thermal_status(heat, target_error, config.reward_config)
    status == thermal_low && return policy.low_action_index
    (status == thermal_high || status == thermal_medium || status == thermal_hard) && return policy.high_action_index
    return zero_action_index()
end
