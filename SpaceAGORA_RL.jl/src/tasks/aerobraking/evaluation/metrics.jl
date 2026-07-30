function episode_metrics(summary::EpisodeSummary; policy_name::AbstractString="")
    return (
        policy = policy_name,
        episode = summary.episode_index,
        worker_id = summary.worker_id,
        seed = summary.seed,
        pass_count = summary.pass_count,
        decision_passes = summary.pass_count - summary.protected_passes,
        protected_passes = summary.protected_passes,
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

function thermal_violation_breakdown(summary::EpisodeSummary,
                                     config::AerobrakingScenarioConfig)
    low = high = medium = hard = 0
    trace_length = min(length(summary.heat_rate_trace), length(summary.apoapsis_trace_m))
    for index in 1:trace_length
        protected = index <= length(summary.protected_trace) &&
                    summary.protected_trace[index]
        protected && continue
        target_error_m = summary.apoapsis_trace_m[index] -
                         config.final_apoapsis_radius_m
        status = thermal_status(
            summary.heat_rate_trace[index],
            target_error_m,
            config.reward_config,
        )
        status == thermal_low && (low += 1)
        status == thermal_high && (high += 1)
        status == thermal_medium && (medium += 1)
        status == thermal_hard && (hard += 1)
    end
    return (
        low = low,
        high = high,
        medium = medium,
        hard = hard,
        total = low + high + medium + hard,
    )
end

function terminal_thermal_violation_type(summary::EpisodeSummary,
                                         config::AerobrakingScenarioConfig)
    config.termination_config.terminal_on_thermal_violation || return nothing
    summary.success && return nothing
    summary.impact && return nothing
    summary.out_of_drag_passage && return nothing
    summary.target_error_m < -config.reward_config.target_tolerance_m && return nothing
    isempty(summary.heat_rate_trace) && return nothing
    last_index = length(summary.heat_rate_trace)
    protected = last_index <= length(summary.protected_trace) &&
                summary.protected_trace[last_index]
    protected && return nothing
    last_index <= length(summary.apoapsis_trace_m) || return nothing
    status = thermal_status(
        summary.heat_rate_trace[last_index],
        summary.apoapsis_trace_m[last_index] - config.final_apoapsis_radius_m,
        config.reward_config,
    )
    return thermal_violation(status) ? status : nothing
end

function episode_thermal_violation_metrics(summary::EpisodeSummary,
                                           config::AerobrakingScenarioConfig)
    counts = thermal_violation_breakdown(summary, config)
    terminal_type = terminal_thermal_violation_type(summary, config)
    terminal_label =
        terminal_type === thermal_low ? "low" :
        terminal_type === thermal_high ? "soft" :
        terminal_type === thermal_medium ? "medium" :
        terminal_type === thermal_hard ? "hard" :
        "none"
    return (
        low_thermal_violations = counts.low,
        soft_thermal_violations = counts.high,
        medium_thermal_violations = counts.medium,
        hard_thermal_violations = counts.hard,
        terminal_thermal_violation_type = terminal_label,
    )
end

function aggregate_thermal_violation_metrics(
    summaries::AbstractVector{<:EpisodeSummary},
    config::AerobrakingScenarioConfig,
)
    isempty(summaries) && return NamedTuple()
    counts = thermal_violation_breakdown.(summaries, Ref(config))
    terminal_types = terminal_thermal_violation_type.(summaries, Ref(config))
    episodes = length(summaries)
    failed = [!summary.success for summary in summaries]

    mean_count(field::Symbol) = sum(getfield(count, field) for count in counts) / episodes
    episode_rate(field::Symbol) =
        count(counts_for_episode -> getfield(counts_for_episode, field) > 0, counts) / episodes
    failed_episode_rate(field::Symbol) =
        count(index -> failed[index] && getfield(counts[index], field) > 0,
              eachindex(summaries)) / episodes
    terminal_rate(status::ThermalStatus) = count(==(status), terminal_types) / episodes

    terminal_low_rate = terminal_rate(thermal_low)
    terminal_high_rate = terminal_rate(thermal_high)
    terminal_medium_rate = terminal_rate(thermal_medium)
    terminal_hard_rate = terminal_rate(thermal_hard)
    return (
        mean_low_thermal_violations = mean_count(:low),
        mean_soft_thermal_violations = mean_count(:high),
        mean_medium_thermal_violations = mean_count(:medium),
        mean_hard_thermal_violations = mean_count(:hard),
        low_violation_episode_rate = episode_rate(:low),
        soft_violation_episode_rate = episode_rate(:high),
        medium_violation_episode_rate = episode_rate(:medium),
        hard_violation_episode_rate = episode_rate(:hard),
        failed_with_low_violation_rate = failed_episode_rate(:low),
        failed_with_soft_violation_rate = failed_episode_rate(:high),
        failed_with_medium_violation_rate = failed_episode_rate(:medium),
        failed_with_hard_violation_rate = failed_episode_rate(:hard),
        terminal_low_violation_rate = terminal_low_rate,
        terminal_soft_violation_rate = terminal_high_rate,
        terminal_medium_violation_rate = terminal_medium_rate,
        terminal_hard_violation_rate = terminal_hard_rate,
        thermal_terminal_failure_rate =
            terminal_low_rate + terminal_high_rate +
            terminal_medium_rate + terminal_hard_rate,
    )
end

function aggregate_metrics(summaries::AbstractVector{<:EpisodeSummary}; policy_name::AbstractString="")
    isempty(summaries) && return NamedTuple()
    successes = 0
    target_error_sum = 0.0
    delta_v_sum = 0.0
    abm_delta_v_sum = 0.0
    pass_count_sum = 0
    protected_passes_sum = 0
    thermal_violations_sum = 0
    solver_failures = 0
    for summary in summaries
        successes += summary.success ? 1 : 0
        target_error_sum += abs(summary.target_error_m)
        delta_v_sum += summary.total_delta_v_mps
        abm_delta_v_sum += summary.abm_delta_v_mps
        pass_count_sum += summary.pass_count
        protected_passes_sum += summary.protected_passes
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
        mean_decision_passes = (pass_count_sum - protected_passes_sum) / episodes,
        mean_protected_passes = protected_passes_sum / episodes,
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
    protected_passes_sum::Int = 0
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
    accumulator.protected_passes_sum += summary.protected_passes
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
        mean_decision_passes =
            (accumulator.pass_count_sum - accumulator.protected_passes_sum) / episodes,
        mean_protected_passes = accumulator.protected_passes_sum / episodes,
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
            protected = summary.protected_trace[i],
        ))
    end
    return rows
end
