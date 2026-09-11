function _safe_stat(values, op::Function)
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return op(vec)
end

@inline function _ci95_bounds(values)
    vec = collect(skipmissing(values))
    n = length(vec)
    if n == 0
        return (missing, missing)
    elseif n == 1
        μ = Float64(vec[1])
        return (μ, μ)
    end
    μ = mean(vec)
    σ = std(vec; corrected=true)
    sem = σ / sqrt(n)
    margin = 1.96 * sem
    return (μ - margin, μ + margin)
end

@inline _ci95_low(values) = _ci95_bounds(values)[1]
@inline _ci95_high(values) = _ci95_bounds(values)[2]

@inline function _sem(values)
    vec = collect(skipmissing(values))
    n = length(vec)
    if n == 0
        return missing
    elseif n == 1
        return 0.0
    end
    return std(vec; corrected=true) / sqrt(n)
end

@inline function _cv_pct(values)
    vec = collect(skipmissing(values))
    n = length(vec)
    if n == 0
        return missing
    elseif n == 1
        return 0.0
    end
    μ = mean(vec)
    abs(μ) <= eps(Float64) && return missing
    return 100.0 * std(vec; corrected=true) / abs(μ)
end

@inline function _sum_nonmissing(values...)
    acc = 0.0
    seen = false
    for v in values
        if !(v isa Missing)
            acc += Float64(v)
            seen = true
        end
    end
    return seen ? acc : missing
end

@inline function _safe_share(num, den)
    if num isa Missing || den isa Missing
        return missing
    end
    n = Float64(num)
    d = Float64(den)
    if !isfinite(n) || !isfinite(d) || d <= 0.0
        return missing
    end
    return n / d
end

@inline function _summary_group_keys(df::DataFrame, base_keys::Vector{Symbol})::Vector{Symbol}
    df_names = Set(Symbol.(names(df)))
    keys = Symbol[]
    for key in base_keys
        if key in df_names
            push!(keys, key)
        end
    end
    optional_keys = (
        :hardware_class,
        :machine_label,
        :host_name,
        :cpu_name,
        :cpu_threads,
        :julia_threads,
        :os,
        :arch
    )
    for key in optional_keys
        if key in df_names && !(key in keys)
            push!(keys, key)
        end
    end
    return keys
end

function summarize_per_orbit_results(orbit_raw_df::DataFrame)::DataFrame
    sweep_multiplier_key = :mission_time_multiplier in names(orbit_raw_df) ? :mission_time_multiplier : :orbit_count
    keys = _summary_group_keys(
        orbit_raw_df,
        [:category, :scenario, :description, sweep_multiplier_key, :orbital_period_s, :dt_max_orbit_s, :outer_threads_safe]
    )
    counts = combine(
        groupby(orbit_raw_df, keys),
        nrow => :samples_total,
        :solve_success => (v -> count(identity, v)) => :samples_success
    )
    counts[!, :samples_failed] = counts.samples_total .- counts.samples_success
    counts[!, :success_rate] = Float64.(counts.samples_success) ./ Float64.(counts.samples_total)

    success_df = orbit_raw_df[orbit_raw_df.solve_success .== true, :]
    summary = counts
    if nrow(success_df) > 0
        success_summary = combine(
            groupby(success_df, keys),
            nrow => :samples,
            :mission_time_s => (v -> _safe_stat(v, mean)) => :mission_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, mean)) => :total_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :total_time_p90_s,
            :total_time_s => (v -> _ci95_low(v)) => :total_time_ci95_low_s,
            :total_time_s => (v -> _ci95_high(v)) => :total_time_ci95_high_s,
            :total_time_s => (v -> _sem(v)) => :total_time_sem_s,
            :total_time_s => (v -> _cv_pct(v)) => :total_time_cv_pct,
            :solve_time_s => (v -> _safe_stat(v, mean)) => :solve_time_mean_s,
            :total_bytes_mb => (v -> _safe_stat(v, mean)) => :total_bytes_mean_mb,
            :sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :sim_seconds_per_wall_second_mean
        )
        summary = leftjoin(counts, success_summary, on=keys, matchmissing=:equal)
    else
        summary[!, :samples] = fill(missing, nrow(summary))
        summary[!, :mission_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_time_p90_s] = fill(missing, nrow(summary))
        summary[!, :total_time_ci95_low_s] = fill(missing, nrow(summary))
        summary[!, :total_time_ci95_high_s] = fill(missing, nrow(summary))
        summary[!, :total_time_sem_s] = fill(missing, nrow(summary))
        summary[!, :total_time_cv_pct] = fill(missing, nrow(summary))
        summary[!, :solve_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_bytes_mean_mb] = fill(missing, nrow(summary))
        summary[!, :sim_seconds_per_wall_second_mean] = fill(missing, nrow(summary))
    end

    if !(:mission_time_multiplier in names(summary))
        summary[!, :mission_time_multiplier] = summary.orbit_count
    end
    summary[!, :time_per_orbit_mean_s] = [
        (ismissing(tt) || multiplier <= 0) ? missing : tt / multiplier
        for (tt, multiplier) in zip(summary.total_time_mean_s, summary.mission_time_multiplier)
    ]
    summary[!, :orbits_per_wall_second_mean] = [
        (ismissing(tt) || tt <= 0.0) ? missing : multiplier / tt
        for (tt, multiplier) in zip(summary.total_time_mean_s, summary.mission_time_multiplier)
    ]
    # Human-facing aliases for mission-time sweep semantics.
    summary[!, :time_per_baseline_period_mean_s] = summary.time_per_orbit_mean_s
    summary[!, :baseline_periods_per_wall_second_mean] = summary.orbits_per_wall_second_mean

    sort!(summary, [:scenario, :mission_time_multiplier])
    return summary
end

function summarize_entry_duration_results(entry_raw_df::DataFrame)::DataFrame
    if nrow(entry_raw_df) == 0
        return DataFrame()
    end

    if !(:entry_run_role in names(entry_raw_df))
        entry_raw_df[!, :entry_run_role] = fill("measured", nrow(entry_raw_df))
    end
    for col in (:entry_atmospheric_interface_count, :entry_passage_duration_s, :entry_wall_time_per_passage_s, :entry_event_time_abs_error_s, :entry_reference_terminal_time_s, :terminal_time_s)
        if !(col in names(entry_raw_df))
            entry_raw_df[!, col] = fill(missing, nrow(entry_raw_df))
        end
    end

    keys = _summary_group_keys(
        entry_raw_df,
        [:category, :scenario, :description, :entry_atmospheric_interface_count, :entry_run_role, :outer_threads_safe]
    )
    counts = combine(
        groupby(entry_raw_df, keys),
        nrow => :samples_total,
        :solve_success => (v -> count(identity, v)) => :samples_success
    )
    counts[!, :samples_failed] = counts.samples_total .- counts.samples_success
    counts[!, :success_rate] = Float64.(counts.samples_success) ./ Float64.(counts.samples_total)

    success_df = entry_raw_df[entry_raw_df.solve_success .== true, :]
    summary = counts
    if nrow(success_df) > 0
        success_summary = combine(
            groupby(success_df, keys),
            nrow => :samples,
            :entry_passage_duration_s => (v -> _safe_stat(v, mean)) => :passage_duration_mean_s,
            :entry_passage_duration_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :passage_duration_p90_s,
            :entry_wall_time_per_passage_s => (v -> _safe_stat(v, mean)) => :wall_time_per_passage_mean_s,
            :entry_wall_time_per_passage_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :wall_time_per_passage_p90_s,
            :entry_event_time_abs_error_s => (v -> _safe_stat(v, mean)) => :event_time_abs_error_mean_s,
            :entry_event_time_abs_error_s => (v -> _safe_stat(v, maximum)) => :event_time_abs_error_max_s,
            :terminal_time_s => (v -> _safe_stat(v, mean)) => :event_time_mean_s,
            :entry_reference_terminal_time_s => (v -> _safe_stat(v, mean)) => :reference_event_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, mean)) => :total_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :total_time_p90_s,
            :sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :sim_seconds_per_wall_second_mean
        )
        summary = leftjoin(counts, success_summary, on=keys, matchmissing=:equal)
    else
        summary[!, :samples] = fill(missing, nrow(summary))
        summary[!, :passage_duration_mean_s] = fill(missing, nrow(summary))
        summary[!, :passage_duration_p90_s] = fill(missing, nrow(summary))
        summary[!, :wall_time_per_passage_mean_s] = fill(missing, nrow(summary))
        summary[!, :wall_time_per_passage_p90_s] = fill(missing, nrow(summary))
        summary[!, :event_time_abs_error_mean_s] = fill(missing, nrow(summary))
        summary[!, :event_time_abs_error_max_s] = fill(missing, nrow(summary))
        summary[!, :event_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :reference_event_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_time_mean_s] = fill(missing, nrow(summary))
        summary[!, :total_time_p90_s] = fill(missing, nrow(summary))
        summary[!, :sim_seconds_per_wall_second_mean] = fill(missing, nrow(summary))
    end

    sort!(summary, [:scenario, :entry_atmospheric_interface_count, :entry_run_role])
    return summary
end

function summarize_results(raw_df::DataFrame)::DataFrame
    keys = _summary_group_keys(
        raw_df,
        [
            :category,
            :scenario,
            :description,
            :satellites,
            :orientation,
            :mission_time_s,
            :outer_route,
            :outer_threads_safe,
            :density_parallel_mode,
            :control_parallel_mode,
            :multibody_parallel_mode,
            :dt_max_orbit_s,
            :dynamic_effectors,
            :control_effectors
        ]
    )
    grouped = groupby(raw_df, keys)
    counts = combine(
        grouped,
        nrow => :samples_total,
        :solve_success => (v -> count(identity, v)) => :samples_success,
        :total_time_s => (v -> begin
            vals = Float64[]
            for x in v
                if !(x isa Missing)
                    push!(vals, Float64(x))
                end
            end
            isempty(vals) ? missing : sum(vals)
        end) => :total_time_all_attempts_s,
        :solver_fallback_count => (v -> begin
            vals = Float64[]
            for x in v
                push!(vals, x isa Missing ? 0.0 : Float64(x))
            end
            isempty(vals) ? missing : mean(vals)
        end) => :fallback_count_mean_all_attempts
    )
    if :attempt in Symbol.(names(raw_df))
        requested_counts = combine(
            grouped,
            :attempt => (v -> begin
                n = 0
                for x in v
                    if !(x isa Missing) && Int(x) == 1
                        n += 1
                    end
                end
                n
            end) => :requested_runs
        )
        counts = leftjoin(counts, requested_counts, on=keys)
    else
        counts[!, :requested_runs] = Int.(counts.samples_total)
    end
    counts[!, :requested_runs] = [
        (v isa Missing || Int(v) <= 0) ? Int(samples_total) : Int(v)
        for (v, samples_total) in zip(counts.requested_runs, counts.samples_total)
    ]
    counts[!, :samples_failed] = counts.samples_total .- counts.samples_success
    counts[!, :success_rate] = Float64.(counts.samples_success) ./ Float64.(counts.samples_total)
    counts[!, :failure_rate] = Float64.(counts.samples_failed) ./ Float64.(counts.samples_total)
    counts[!, :retries_total] = max.(0, counts.samples_total .- counts.requested_runs)
    counts[!, :retry_count_mean] = [
        requested <= 0 ? missing : (Float64(retries) / Float64(requested))
        for (retries, requested) in zip(counts.retries_total, counts.requested_runs)
    ]
    counts[!, :penalized_expected_wall_time_s] = [
        (total_time isa Missing || requested <= 0) ? missing : (Float64(total_time) / Float64(requested))
        for (total_time, requested) in zip(counts.total_time_all_attempts_s, counts.requested_runs)
    ]

    success_df = raw_df[raw_df.solve_success .== true, :]
    metric_cols = [
        :samples,
        :copy_time_mean_s,
        :copy_bytes_mean_mb,
        :solve_time_mean_s,
        :solve_bytes_mean_mb,
        :total_time_mean_s,
        :copy_compile_time_mean_s,
        :solve_compile_time_mean_s,
        :compile_time_mean_s,
        :copy_gc_time_mean_s,
        :solve_gc_time_mean_s,
        :gc_time_mean_s,
        :setup_share,
        :copy_time_share,
        :solve_share,
        :compile_share,
        :gc_share,
        :compile_gc_share,
        :copy_bytes_share,
        :copy_alloc_share,
        :total_time_median_s,
        :total_time_std_s,
        :total_time_min_s,
        :total_time_max_s,
        :total_time_p90_s,
        :total_time_ci95_low_s,
        :total_time_ci95_high_s,
        :total_time_sem_s,
        :total_time_cv_pct,
        :total_bytes_mean_mb,
        :copy_alloc_mean,
        :solve_alloc_mean,
        :saved_points_mean,
        :accepted_steps_mean,
        :rejected_steps_mean,
        :solver_modes,
        :solver_sequences,
        :solver_fallback_any,
        :solver_fallback_count_mean,
        :solver_fallback_triggers,
        :parallel_threads_used_any,
        :policy_decisions_mean,
        :policy_threads_enabled_mean,
        :policy_density_threads_enabled_mean,
        :policy_control_threads_enabled_mean,
        :policy_multibody_threads_enabled_mean,
        :policy_other_threads_enabled_mean,
        :sim_seconds_per_wall_second_mean,
        :satellite_sim_seconds_per_wall_second_mean,
        :nbody_spkpos_runtime_calls_mean,
        :nbody_spkpos_cache_build_calls_mean,
        :nbody_spkpos_total_calls_mean,
        :srp_spkpos_runtime_calls_mean,
        :srp_spkpos_cache_build_calls_mean,
        :srp_spkpos_total_calls_mean,
        :planet_pxform_runtime_calls_mean,
        :planet_pxform_cache_build_calls_mean,
        :planet_pxform_total_calls_mean,
        :spice_calls_runtime_mean,
        :spice_calls_cache_build_mean,
        :spice_calls_total_mean,
        :spice_calls_per_wall_second_mean,
        :spice_calls_per_sim_second_mean
    ]

    summary = counts
    if nrow(success_df) > 0
        success_summary = combine(
            groupby(success_df, keys),
            nrow => :samples,
            :copy_time_s => (v -> _safe_stat(v, mean)) => :copy_time_mean_s,
            :copy_bytes_mb => (v -> _safe_stat(v, mean)) => :copy_bytes_mean_mb,
            :solve_time_s => (v -> _safe_stat(v, mean)) => :solve_time_mean_s,
            :solve_bytes_mb => (v -> _safe_stat(v, mean)) => :solve_bytes_mean_mb,
            :total_time_s => (v -> _safe_stat(v, mean)) => :total_time_mean_s,
            :copy_compile_time_s => (v -> _safe_stat(v, mean)) => :copy_compile_time_mean_s,
            :solve_compile_time_s => (v -> _safe_stat(v, mean)) => :solve_compile_time_mean_s,
            :copy_gctime_s => (v -> _safe_stat(v, mean)) => :copy_gc_time_mean_s,
            :solve_gctime_s => (v -> _safe_stat(v, mean)) => :solve_gc_time_mean_s,
            :total_time_s => (v -> _safe_stat(v, median)) => :total_time_median_s,
            :total_time_s => (v -> _safe_stat(v, x -> std(x; corrected=false))) => :total_time_std_s,
            :total_time_s => (v -> _safe_stat(v, minimum)) => :total_time_min_s,
            :total_time_s => (v -> _safe_stat(v, maximum)) => :total_time_max_s,
            :total_time_s => (v -> _safe_stat(v, x -> quantile(x, 0.9))) => :total_time_p90_s,
            :total_time_s => (v -> _ci95_low(v)) => :total_time_ci95_low_s,
            :total_time_s => (v -> _ci95_high(v)) => :total_time_ci95_high_s,
            :total_time_s => (v -> _sem(v)) => :total_time_sem_s,
            :total_time_s => (v -> _cv_pct(v)) => :total_time_cv_pct,
            :total_bytes_mb => (v -> _safe_stat(v, mean)) => :total_bytes_mean_mb,
            :copy_alloc_calls => (v -> _safe_stat(v, mean)) => :copy_alloc_mean,
            :solve_alloc_calls => (v -> _safe_stat(v, mean)) => :solve_alloc_mean,
            :saved_points => (v -> _safe_stat(v, mean)) => :saved_points_mean,
            :accepted_steps => (v -> _safe_stat(v, mean)) => :accepted_steps_mean,
            :rejected_steps => (v -> _safe_stat(v, mean)) => :rejected_steps_mean,
            :solver_mode => (v -> _safe_unique_join(v)) => :solver_modes,
            :solver_sequence => (v -> _safe_unique_join(v)) => :solver_sequences,
            :solver_fallback_used => (v -> any(skipmissing(v))) => :solver_fallback_any,
            :solver_fallback_count => (v -> _safe_stat(v, mean)) => :solver_fallback_count_mean,
            :solver_fallback_trigger => (v -> _safe_unique_join(v; delimiter="|")) => :solver_fallback_triggers,
            :policy_threads_enabled_total => (v -> begin
                vals = collect(skipmissing(v))
                isempty(vals) ? missing : any(x -> x > 0, vals)
            end) => :parallel_threads_used_any,
            :policy_decisions_total => (v -> _safe_stat(v, mean)) => :policy_decisions_mean,
            :policy_threads_enabled_total => (v -> _safe_stat(v, mean)) => :policy_threads_enabled_mean,
            :policy_density_threads_enabled => (v -> _safe_stat(v, mean)) => :policy_density_threads_enabled_mean,
            :policy_control_threads_enabled => (v -> _safe_stat(v, mean)) => :policy_control_threads_enabled_mean,
            :policy_multibody_threads_enabled => (v -> _safe_stat(v, mean)) => :policy_multibody_threads_enabled_mean,
            :policy_other_threads_enabled => (v -> _safe_stat(v, mean)) => :policy_other_threads_enabled_mean,
            :sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :sim_seconds_per_wall_second_mean,
            :satellite_sim_seconds_per_wall_second => (v -> _safe_stat(v, mean)) => :satellite_sim_seconds_per_wall_second_mean,
            :nbody_spkpos_runtime_calls => (v -> _safe_stat(v, mean)) => :nbody_spkpos_runtime_calls_mean,
            :nbody_spkpos_cache_build_calls => (v -> _safe_stat(v, mean)) => :nbody_spkpos_cache_build_calls_mean,
            :nbody_spkpos_total_calls => (v -> _safe_stat(v, mean)) => :nbody_spkpos_total_calls_mean,
            :srp_spkpos_runtime_calls => (v -> _safe_stat(v, mean)) => :srp_spkpos_runtime_calls_mean,
            :srp_spkpos_cache_build_calls => (v -> _safe_stat(v, mean)) => :srp_spkpos_cache_build_calls_mean,
            :srp_spkpos_total_calls => (v -> _safe_stat(v, mean)) => :srp_spkpos_total_calls_mean,
            :planet_pxform_runtime_calls => (v -> _safe_stat(v, mean)) => :planet_pxform_runtime_calls_mean,
            :planet_pxform_cache_build_calls => (v -> _safe_stat(v, mean)) => :planet_pxform_cache_build_calls_mean,
            :planet_pxform_total_calls => (v -> _safe_stat(v, mean)) => :planet_pxform_total_calls_mean
        )
        summary = leftjoin(counts, success_summary, on=keys)
    else
        for col in metric_cols
            summary[!, col] = fill(missing, nrow(summary))
        end
    end

    summary[!, :compile_time_mean_s] = [
        _sum_nonmissing(row.copy_compile_time_mean_s, row.solve_compile_time_mean_s)
        for row in eachrow(summary)
    ]
    summary[!, :gc_time_mean_s] = [
        _sum_nonmissing(row.copy_gc_time_mean_s, row.solve_gc_time_mean_s)
        for row in eachrow(summary)
    ]
    summary[!, :setup_share] = [
        _safe_share(row.copy_time_mean_s, row.total_time_mean_s)
        for row in eachrow(summary)
    ]
    summary[!, :copy_time_share] = copy(summary.setup_share)
    summary[!, :solve_share] = [
        _safe_share(row.solve_time_mean_s, row.total_time_mean_s)
        for row in eachrow(summary)
    ]
    summary[!, :compile_share] = [
        _safe_share(row.compile_time_mean_s, row.total_time_mean_s)
        for row in eachrow(summary)
    ]
    summary[!, :gc_share] = [
        _safe_share(row.gc_time_mean_s, row.total_time_mean_s)
        for row in eachrow(summary)
    ]
    summary[!, :compile_gc_share] = [
        _safe_share(_sum_nonmissing(row.compile_time_mean_s, row.gc_time_mean_s), row.total_time_mean_s)
        for row in eachrow(summary)
    ]
    summary[!, :copy_bytes_share] = [
        _safe_share(row.copy_bytes_mean_mb, row.total_bytes_mean_mb)
        for row in eachrow(summary)
    ]
    summary[!, :copy_alloc_share] = [
        _safe_share(row.copy_alloc_mean, _sum_nonmissing(row.copy_alloc_mean, row.solve_alloc_mean))
        for row in eachrow(summary)
    ]

    summary[!, :spice_calls_runtime_mean] = [
        _sum_nonmissing(row.nbody_spkpos_runtime_calls_mean, row.srp_spkpos_runtime_calls_mean, row.planet_pxform_runtime_calls_mean)
        for row in eachrow(summary)
    ]
    summary[!, :spice_calls_cache_build_mean] = [
        _sum_nonmissing(row.nbody_spkpos_cache_build_calls_mean, row.srp_spkpos_cache_build_calls_mean, row.planet_pxform_cache_build_calls_mean)
        for row in eachrow(summary)
    ]
    summary[!, :spice_calls_total_mean] = [
        _sum_nonmissing(row.nbody_spkpos_total_calls_mean, row.srp_spkpos_total_calls_mean, row.planet_pxform_total_calls_mean)
        for row in eachrow(summary)
    ]
    summary[!, :spice_calls_per_wall_second_mean] = [
        (row.spice_calls_total_mean isa Missing || row.total_time_mean_s isa Missing || row.total_time_mean_s <= 0.0) ? missing :
        (Float64(row.spice_calls_total_mean) / Float64(row.total_time_mean_s))
        for row in eachrow(summary)
    ]
    summary[!, :spice_calls_per_sim_second_mean] = [
        (row.spice_calls_total_mean isa Missing || row.mission_time_s isa Missing || row.mission_time_s <= 0.0) ? missing :
        (Float64(row.spice_calls_total_mean) / Float64(row.mission_time_s))
        for row in eachrow(summary)
    ]

    baseline_idx = findfirst(==(PERF_BASELINE_SCENARIO), summary.scenario)
    if baseline_idx === nothing || ismissing(summary.total_time_mean_s[baseline_idx]) || summary.total_time_mean_s[baseline_idx] <= 0.0
        summary[!, :relative_to_baseline] = fill(missing, nrow(summary))
        summary[!, :speedup_vs_baseline] = fill(missing, nrow(summary))
    else
        baseline = summary.total_time_mean_s[baseline_idx]
        summary[!, :relative_to_baseline] = [
            ismissing(row_total) ? missing : row_total / baseline
            for row_total in summary.total_time_mean_s
        ]
        summary[!, :speedup_vs_baseline] = [
            ismissing(row_total) || row_total <= 0.0 ? missing : baseline / row_total
            for row_total in summary.total_time_mean_s
        ]
    end

    summary[!, :_sort_key] = [ismissing(x) ? -Inf : x for x in summary.total_time_mean_s]
    sort!(summary, :_sort_key, rev=true)
    select!(summary, Not(:_sort_key))
    return summary
end

function summarize_density_backend_breakdown(raw_df::DataFrame)::DataFrame
    expected_buckets = [
        "gram_point_to_point",
        "gram_surrogate",
        "gram_static_grid_or_cached_surrogate",
        "non_gram"
    ]
    if nrow(raw_df) == 0 || !(:density_backend_bucket in names(raw_df))
        return DataFrame([
            (
                density_backend_bucket=bucket,
                covered=false,
                density_families="",
                outer_routes="",
                samples_total=0,
                samples_success=0,
                samples_failed=0,
                success_rate_pct=missing,
                total_time_mean_s=missing,
                total_time_p90_s=missing,
                sim_seconds_per_wall_second_mean=missing,
                recommended_route=(bucket == "gram_point_to_point" ? "process" : "threads_or_auto")
            ) for bucket in expected_buckets
        ])
    end

    rows = NamedTuple[]
    for bucket in expected_buckets
        bucket_df = raw_df[raw_df.density_backend_bucket .== bucket, :]
        samples_total = nrow(bucket_df)
        samples_success = samples_total == 0 ? 0 : count(identity, bucket_df.solve_success)
        samples_failed = samples_total - samples_success
        success_rate_pct = samples_total > 0 ? (100.0 * samples_success / samples_total) : missing
        success_df = samples_total == 0 ? bucket_df : bucket_df[bucket_df.solve_success .== true, :]
        total_time_mean_s = nrow(success_df) > 0 ? _safe_stat(success_df.total_time_s, mean) : missing
        total_time_p90_s = nrow(success_df) > 0 ? _safe_stat(success_df.total_time_s, x -> quantile(x, 0.9)) : missing
        throughput_mean = nrow(success_df) > 0 ? _safe_stat(success_df.sim_seconds_per_wall_second, mean) : missing
        push!(rows, (
            density_backend_bucket=bucket,
            covered=samples_total > 0,
            density_families=_safe_unique_join(bucket_df.density_family),
            outer_routes=_safe_unique_join(bucket_df.outer_route),
            samples_total=samples_total,
            samples_success=samples_success,
            samples_failed=samples_failed,
            success_rate_pct=success_rate_pct,
            total_time_mean_s=total_time_mean_s,
            total_time_p90_s=total_time_p90_s,
            sim_seconds_per_wall_second_mean=throughput_mean,
            recommended_route=(bucket == "gram_point_to_point" ? "process" : "threads_or_auto")
        ))
    end
    return DataFrame(rows)
end

@inline function _case_with_solver(
    case::BenchmarkCase;
    solver_mode_override::Union{Nothing, String}=nothing,
    split_imex_solver_override::Union{Nothing, String}=nothing
)::BenchmarkCase
    return BenchmarkCase(
        name=case.name,
        category=case.category,
        description=case.description,
        args_template=case.args_template,
        run_in_quick=case.run_in_quick,
        solver_mode_override=solver_mode_override,
        split_imex_solver_override=split_imex_solver_override,
        entry_target_count_override=case.entry_target_count_override,
        env_overrides=copy(case.env_overrides)
    )
end

@inline function _solution_state_at(sol, t::Float64)
    try
        return sol(t)
    catch
        idx = searchsortedlast(sol.t, t)
        idx = clamp(idx, 1, length(sol.u))
        return sol.u[idx]
    end
end

@inline function _relative_vector_delta(ref_vec, cmp_vec; floor::Float64=1e-12)::Float64
    ref_norm = norm(ref_vec)
    delta_norm = norm(ref_vec - cmp_vec)
    return ref_norm > floor ? delta_norm / ref_norm : delta_norm
end

@inline function _quaternion_angle_delta_rad(q_ref, q_cmp)::Float64
    q_ref_norm = norm(q_ref)
    q_cmp_norm = norm(q_cmp)
    if !(isfinite(q_ref_norm) && isfinite(q_cmp_norm)) || q_ref_norm <= eps(Float64) || q_cmp_norm <= eps(Float64)
        return Inf
    end
    c = abs(dot(q_ref / q_ref_norm, q_cmp / q_cmp_norm))
    c = clamp(c, -1.0, 1.0)
    return 2.0 * acos(c)
end
