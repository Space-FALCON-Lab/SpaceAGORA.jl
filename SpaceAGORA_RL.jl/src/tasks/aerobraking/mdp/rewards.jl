Base.@kwdef struct RewardConfig
    target_tolerance_m::Float64 = 10e3
    near_target_band_m::Float64 = 100e3
    heat_low_w_cm2::Float64 = 0.05
    heat_high_w_cm2::Float64 = 0.25
    heat_medium_w_cm2::Float64 = 0.30
    heat_hard_w_cm2::Float64 = 0.45
    impact_penalty::Float64 = -6.0
    out_of_passage_penalty::Float64 = -6.0
    undershoot_penalty::Float64 = -4.0
    low_thermal_penalty::Float64 = -4.0
    medium_thermal_penalty::Float64 = -5.0
    hard_thermal_penalty::Float64 = -6.0
    success_reward_base::Float64 = 10.0
    success_distance_scale::Float64 = 2 / 5
    gaussian_scale::Float64 = 0.025
    heat_center::Float64 = 0.5
    heat_sigma::Float64 = 0.1
    step_penalty::Float64 = -0.1
    action_penalty_scale::Float64 = 0.1
    low_heat_action_cap_mps::Float64 = 1.5
end

@enum ThermalStatus thermal_low thermal_nominal thermal_high thermal_medium thermal_hard

function thermal_status(max_heat_rate_w_cm2::Real, target_error_m::Real, config::RewardConfig)
    heat = Float64(max_heat_rate_w_cm2)
    far_from_target = abs(Float64(target_error_m)) > config.near_target_band_m
    if heat > config.heat_hard_w_cm2
        return thermal_hard
    elseif heat > config.heat_medium_w_cm2
        return thermal_medium
    elseif heat > config.heat_high_w_cm2
        return thermal_high
    elseif far_from_target && heat < config.heat_low_w_cm2
        return thermal_low
    else
        return thermal_nominal
    end
end

thermal_violation(status::ThermalStatus) = status != thermal_nominal

function thermal_penalty(status::ThermalStatus, config::RewardConfig)
    status == thermal_hard && return config.hard_thermal_penalty
    status == thermal_medium && return config.medium_thermal_penalty
    return config.low_thermal_penalty
end

function paper_reward(obs::PaperObservation, scenario_config, action::AerobrakingAction,
                      flags, config::RewardConfig=scenario_config.reward_config)
    ra = obs.apoapsis_radius_m
    target = scenario_config.final_apoapsis_radius_m
    heat = obs.max_heat_rate_w_cm2
    reward = 0.0

    if flags.terminated
        if flags.impact
            reward += config.impact_penalty
        elseif flags.out_of_drag_passage
            reward += config.out_of_passage_penalty
        elseif flags.success
            dist_km_floor = floor(abs(ra - target) / 1000)
            reward += -dist_km_floor * config.success_distance_scale + config.success_reward_base
        elseif flags.target_undershoot
            reward += config.undershoot_penalty
        end
    end

    status = thermal_status(heat, ra - target, config)
    if thermal_violation(status)
        reward += thermal_penalty(status, config)
    else
        heat_cap = (heat - config.heat_low_w_cm2) / (config.heat_high_w_cm2 - config.heat_low_w_cm2)
        delta_v_cap = abs(action.delta_v_mps)
        if heat_cap < 0
            heat_cap = 0.0
            delta_v_cap = abs(action.magnitude_mps - config.low_heat_action_cap_mps)
        end
        distance_cap = (ra - target) / (scenario_config.r_norm_m - target)
        if distance_cap >= 0
            thermal_shape = config.gaussian_scale / (config.heat_sigma * sqrt(2pi)) *
                            exp(-0.5 * ((heat_cap - config.heat_center) / config.heat_sigma)^2)
            reward += (1 - distance_cap^0.4) * thermal_shape +
                      (config.step_penalty - config.action_penalty_scale * delta_v_cap)
        end
    end

    return reward
end
