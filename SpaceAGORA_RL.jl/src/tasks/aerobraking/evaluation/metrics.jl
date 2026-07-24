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

function aggregate_metrics(summaries::AbstractVector{<:EpisodeSummary}; policy_name::AbstractString="")
    isempty(summaries) && return NamedTuple()
    successes = 0
    target_error_sum = 0.0
    delta_v_sum = 0.0
    abm_delta_v_sum = 0.0
    pass_count_sum = 0
    thermal_violations_sum = 0
    solver_failures = 0
    for summary in summaries
        successes += summary.success ? 1 : 0
        target_error_sum += abs(summary.target_error_m)
        delta_v_sum += summary.total_delta_v_mps
        abm_delta_v_sum += summary.abm_delta_v_mps
        pass_count_sum += summary.pass_count
        thermal_violations_sum += summary.thermal_violations
        solver_failures += summary.solver_failures
    end
    episodes = length(summaries)
    return (
        policy = policy_name,
        episodes = episodes,
        success_rate = successes / episodes,
        mean_target_error_km = target_error_sum / episodes / 1000,
        mean_delta_v_mps = delta_v_sum / episodes,
        mean_abm_delta_v_mps = abm_delta_v_sum / episodes,
        mean_pass_count = pass_count_sum / episodes,
        mean_thermal_violations = thermal_violations_sum / episodes,
        solver_failures = solver_failures,
    )
end

Base.@kwdef mutable struct EpisodeAggregateAccumulator
    episodes::Int = 0
    successes::Int = 0
    target_error_sum::Float64 = 0.0
    delta_v_sum::Float64 = 0.0
    abm_delta_v_sum::Float64 = 0.0
    pass_count_sum::Int = 0
    thermal_violations_sum::Int = 0
    solver_failures::Int = 0
end

function accumulate_episode!(accumulator::EpisodeAggregateAccumulator, summary::EpisodeSummary)
    accumulator.episodes += 1
    accumulator.successes += summary.success ? 1 : 0
    accumulator.target_error_sum += abs(summary.target_error_m)
    accumulator.delta_v_sum += summary.total_delta_v_mps
    accumulator.abm_delta_v_sum += summary.abm_delta_v_mps
    accumulator.pass_count_sum += summary.pass_count
    accumulator.thermal_violations_sum += summary.thermal_violations
    accumulator.solver_failures += summary.solver_failures
    return accumulator
end

function aggregate_metrics(accumulator::EpisodeAggregateAccumulator; policy_name::AbstractString="")
    accumulator.episodes == 0 && return NamedTuple()
    episodes = accumulator.episodes
    return (
        policy = policy_name,
        episodes = episodes,
        success_rate = accumulator.successes / episodes,
        mean_target_error_km = accumulator.target_error_sum / episodes / 1000,
        mean_delta_v_mps = accumulator.delta_v_sum / episodes,
        mean_abm_delta_v_mps = accumulator.abm_delta_v_sum / episodes,
        mean_pass_count = accumulator.pass_count_sum / episodes,
        mean_thermal_violations = accumulator.thermal_violations_sum / episodes,
        solver_failures = accumulator.solver_failures,
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
