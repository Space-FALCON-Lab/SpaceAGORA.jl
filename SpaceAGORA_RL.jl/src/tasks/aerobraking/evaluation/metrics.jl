function episode_metrics(summary::EpisodeSummary; policy_name::AbstractString="")
    return (
        policy = policy_name,
        episode = summary.episode_index,
        worker_id = summary.worker_id,
        seed = summary.seed,
        pass_count = summary.pass_count,
        episode_reward = summary.episode_reward,
        success = summary.success,
        impact = summary.impact,
        out_of_drag_passage = summary.out_of_drag_passage,
        thermal_violations = summary.thermal_violations,
        target_error_km = summary.target_error_m / 1000,
        mission_duration_days = summary.mission_duration_days,
        total_delta_v_mps = summary.total_delta_v_mps,
        abm_delta_v_mps = summary.abm_delta_v_mps,
        apoapsis_correction_delta_v_mps = summary.apoapsis_correction_delta_v_mps,
        periapsis_raise_delta_v_mps = summary.periapsis_raise_delta_v_mps,
        total_mission_delta_v_mps = summary.total_mission_delta_v_mps,
        maneuver_count = summary.maneuver_count,
        peak_heat_rate_w_cm2 = isempty(summary.heat_rate_trace) ? NaN : maximum(summary.heat_rate_trace),
        minimum_periapsis_km = isempty(summary.periapsis_trace_m) ? NaN : minimum(summary.periapsis_trace_m) / 1000,
        solver_failures = summary.solver_failures,
    )
end

function aggregate_metrics(summaries::Vector{EpisodeSummary}; policy_name::AbstractString="")
    isempty(summaries) && return NamedTuple()
    success_rate = count(s -> s.success, summaries) / length(summaries)
    return (
        policy = policy_name,
        episodes = length(summaries),
        success_rate = success_rate,
        mean_target_error_km = mean(abs.(getfield.(summaries, :target_error_m))) / 1000,
        mean_delta_v_mps = mean(getfield.(summaries, :total_delta_v_mps)),
        mean_abm_delta_v_mps = mean(getfield.(summaries, :abm_delta_v_mps)),
        mean_pass_count = mean(getfield.(summaries, :pass_count)),
        mean_thermal_violations = mean(getfield.(summaries, :thermal_violations)),
        solver_failures = sum(getfield.(summaries, :solver_failures)),
    )
end

function pass_log_rows(summary::EpisodeSummary; policy_name::AbstractString="")
    rows = NamedTuple[]
    for i in 1:summary.pass_count
        push!(rows, (
            policy = policy_name,
            episode = summary.episode_index,
            pass = i,
            heat_rate_w_cm2 = summary.heat_rate_trace[i],
            action_delta_v_mps = summary.action_trace[i],
            apoapsis_radius_km = summary.apoapsis_trace_m[i] / 1000,
            periapsis_altitude_km = summary.periapsis_trace_m[i] / 1000,
            argument_of_periapsis_rad = summary.omega_trace_rad[i],
            raan_rad = summary.raan_trace_rad[i],
            inclination_rad = summary.inclination_trace_rad[i],
            reward = summary.reward_trace[i],
        ))
    end
    return rows
end
