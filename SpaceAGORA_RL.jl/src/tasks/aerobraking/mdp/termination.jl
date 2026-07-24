Base.@kwdef struct TerminationConfig
    impact_periapsis_altitude_m::Float64 = 85e3
    out_of_passage_periapsis_altitude_m::Float64 = 145e3
    max_passes::Int = 1000
    terminal_on_thermal_violation::Bool = true
end

struct TerminationFlags
    success::Bool
    target_undershoot::Bool
    impact::Bool
    out_of_drag_passage::Bool
    thermal_violation::Bool
    terminated::Bool
    truncated::Bool
end

function classify_termination(obs::PaperObservation, scenario_config,
                              reward_config::RewardConfig=scenario_config.reward_config,
                              termination_config::TerminationConfig=scenario_config.termination_config;
                              training::Bool=true,
                              pass_count::Integer=0)
    ra = obs.apoapsis_radius_m
    hp = obs.periapsis_altitude_m
    target = scenario_config.final_apoapsis_radius_m
    success = abs(ra - target) <= reward_config.target_tolerance_m
    target_undershoot = (ra - target) < -reward_config.target_tolerance_m
    impact = hp < termination_config.impact_periapsis_altitude_m
    out_of_drag_passage = hp > termination_config.out_of_passage_periapsis_altitude_m
    status = thermal_status(obs.max_heat_rate_w_cm2, ra - target, reward_config)
    thermal_bad = thermal_violation(status)
    thermal_terminal = termination_config.terminal_on_thermal_violation && thermal_bad
    terminated = success || target_undershoot || impact || out_of_drag_passage || thermal_terminal
    truncated = pass_count >= termination_config.max_passes
    return TerminationFlags(success, target_undershoot, impact, out_of_drag_passage,
                            thermal_bad, terminated, truncated)
end
