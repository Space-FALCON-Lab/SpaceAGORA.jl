function transition_from_step(previous_obs::Vector{Float32}, action_index::Int,
                              result::AerobrakingStepResult, info_index::Int)
    return Transition(previous_obs, action_index, Float32(result.reward),
                      result.normalized_observation, result.flags.terminated,
                      result.flags.truncated, info_index)
end

function update_episode_summary(summary::EpisodeSummary, result::AerobrakingStepResult)
    metrics = result.metrics
    summary.pass_count = metrics.pass_index
    summary.episode_reward += result.reward
    summary.success = result.flags.success
    summary.impact = result.flags.impact
    summary.out_of_drag_passage = result.flags.out_of_drag_passage
    summary.thermal_violations += result.flags.thermal_violation ? 1 : 0
    summary.target_error_m = NaN
    summary.mission_duration_days = metrics.mission_elapsed_s / 86400
    summary.total_delta_v_mps = metrics.total_delta_v_mps
    summary.abm_delta_v_mps += result.action.magnitude_mps
    summary.apoapsis_correction_delta_v_mps = 0.0
    summary.periapsis_raise_delta_v_mps = 0.0
    summary.total_mission_delta_v_mps = metrics.total_delta_v_mps
    summary.maneuver_count = metrics.maneuver_count
    summary.solver_failures += metrics.solver_converged ? 0 : 1
    push!(summary.heat_rate_trace, metrics.max_heat_rate_w_cm2)
    push!(summary.action_trace, result.action.delta_v_mps)
    push!(summary.apoapsis_trace_m, metrics.apoapsis_radius_m)
    push!(summary.periapsis_trace_m, metrics.periapsis_altitude_m)
    push!(summary.omega_trace_rad, metrics.argument_of_periapsis_rad)
    push!(summary.raan_trace_rad, metrics.raan_rad)
    push!(summary.inclination_trace_rad, metrics.inclination_rad)
    push!(summary.reward_trace, result.reward)
    return summary
end

function finalize_episode_summary(summary::EpisodeSummary, config)
    target_error = isempty(summary.apoapsis_trace_m) ? NaN :
                   last(summary.apoapsis_trace_m) - config.final_apoapsis_radius_m
    return EpisodeSummary(
        episode_index = summary.episode_index,
        worker_id = summary.worker_id,
        seed = summary.seed,
        pass_count = summary.pass_count,
        episode_reward = summary.episode_reward,
        success = summary.success,
        impact = summary.impact,
        out_of_drag_passage = summary.out_of_drag_passage,
        thermal_violations = summary.thermal_violations,
        target_error_m = target_error,
        mission_duration_days = summary.mission_duration_days,
        total_delta_v_mps = summary.total_delta_v_mps,
        abm_delta_v_mps = summary.abm_delta_v_mps,
        apoapsis_correction_delta_v_mps = summary.apoapsis_correction_delta_v_mps,
        periapsis_raise_delta_v_mps = summary.periapsis_raise_delta_v_mps,
        total_mission_delta_v_mps = summary.total_mission_delta_v_mps,
        maneuver_count = summary.maneuver_count,
        solver_failures = summary.solver_failures,
        heat_rate_trace = summary.heat_rate_trace,
        action_trace = summary.action_trace,
        apoapsis_trace_m = summary.apoapsis_trace_m,
        periapsis_trace_m = summary.periapsis_trace_m,
        omega_trace_rad = summary.omega_trace_rad,
        raan_trace_rad = summary.raan_trace_rad,
        inclination_trace_rad = summary.inclination_trace_rad,
        reward_trace = summary.reward_trace,
    )
end
