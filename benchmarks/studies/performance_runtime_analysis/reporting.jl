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
        :solve_time_mean_s,
        :total_time_mean_s,
        :copy_compile_time_mean_s,
        :solve_compile_time_mean_s,
        :compile_time_mean_s,
        :copy_gc_time_mean_s,
        :solve_gc_time_mean_s,
        :gc_time_mean_s,
        :setup_share,
        :solve_share,
        :compile_share,
        :gc_share,
        :compile_gc_share,
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
            :solve_time_s => (v -> _safe_stat(v, mean)) => :solve_time_mean_s,
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

function _trajectory_delta_metrics(
    sol_ref,
    sol_cmp,
    n_sats::Int,
    orientation::Bool;
    n_samples::Int
)
    t_ref_start = Float64(first(sol_ref.t))
    t_ref_end = Float64(last(sol_ref.t))
    t_cmp_start = Float64(first(sol_cmp.t))
    t_cmp_end = Float64(last(sol_cmp.t))
    t_start = max(t_ref_start, t_cmp_start)
    t_end = min(t_ref_end, t_cmp_end)
    if !(isfinite(t_start) && isfinite(t_end)) || t_end < t_start
        throw(ArgumentError("Cannot compare trajectories: incompatible time spans [$(t_ref_start), $(t_ref_end)] vs [$(t_cmp_start), $(t_cmp_end)]."))
    end
    sample_count = max(2, n_samples)
    sample_times = collect(range(t_start, t_end; length=sample_count))

    pos_rel_max = 0.0
    vel_rel_max = 0.0
    q_angle_max_rad = 0.0
    omega_rel_max = 0.0
    pos_rel_sum2 = 0.0
    vel_rel_sum2 = 0.0
    q_angle_sum2 = 0.0
    omega_rel_sum2 = 0.0
    sample_count_total = 0
    orientation_count_total = 0

    for t in sample_times
        u_ref = _solution_state_at(sol_ref, t)
        u_cmp = _solution_state_at(sol_cmp, t)
        for sat_idx in 1:n_sats
            sc_ref = u_ref.sc[sat_idx]
            sc_cmp = u_cmp.sc[sat_idx]
            pos_rel = _relative_vector_delta(sc_ref.pos, sc_cmp.pos)
            vel_rel = _relative_vector_delta(sc_ref.vel, sc_cmp.vel)
            pos_rel_max = max(pos_rel_max, pos_rel)
            vel_rel_max = max(vel_rel_max, vel_rel)
            pos_rel_sum2 += pos_rel^2
            vel_rel_sum2 += vel_rel^2
            sample_count_total += 1
            if orientation
                q_angle = _quaternion_angle_delta_rad(sc_ref.q, sc_cmp.q)
                omega_rel = _relative_vector_delta(sc_ref.ω, sc_cmp.ω)
                q_angle_max_rad = max(q_angle_max_rad, q_angle)
                omega_rel_max = max(omega_rel_max, omega_rel)
                q_angle_sum2 += q_angle^2
                omega_rel_sum2 += omega_rel^2
                orientation_count_total += 1
            end
        end
    end

    return (
        t_start=t_start,
        t_end=t_end,
        sample_count=sample_count,
        compared_state_count=sample_count_total,
        pos_rel_max=pos_rel_max,
        pos_rel_rms=sample_count_total > 0 ? sqrt(pos_rel_sum2 / sample_count_total) : missing,
        vel_rel_max=vel_rel_max,
        vel_rel_rms=sample_count_total > 0 ? sqrt(vel_rel_sum2 / sample_count_total) : missing,
        q_angle_max_rad=orientation ? q_angle_max_rad : missing,
        q_angle_rms_rad=orientation && orientation_count_total > 0 ? sqrt(q_angle_sum2 / orientation_count_total) : missing,
        omega_rel_max=orientation ? omega_rel_max : missing,
        omega_rel_rms=orientation && orientation_count_total > 0 ? sqrt(omega_rel_sum2 / orientation_count_total) : missing
    )
end

function _run_split_gate_solution(
    case::BenchmarkCase,
    profile_name::String;
    solver_mode::String,
    split_solver::Union{Nothing, String}=nothing
)
    GC.gc()
    args_run = deepcopy(case.args_template)
    started_ns = time_ns()
    try
        run_once = () -> _run_perf_simulation(
            args_run;
            return_solution=true,
            return_solver_metadata=true,
            profile_name=profile_name,
            solver_mode_override=solver_mode,
            split_imex_solver_override=split_solver,
            entry_target_count_override=case.entry_target_count_override
        )
        result = isempty(case.env_overrides) ? run_once() : withenv(case.env_overrides...) do
            run_once()
        end
        elapsed_s = (time_ns() - started_ns) / 1e9
        sol = result.solution
        solver_trace = result.solver_trace
        solver_sequence = isempty(solver_trace) ? missing : join([meta.solver for meta in solver_trace], "->")
        return (
            ok=true,
            elapsed_s=elapsed_s,
            success=_solve_success_for_case(sol, case),
            retcode=string(sol.retcode),
            solver_mode=result.solver_mode,
            solver_sequence=solver_sequence,
            solution=sol,
            error_text=missing
        )
    catch err
        if err isa InterruptException
            rethrow()
        end
        elapsed_s = (time_ns() - started_ns) / 1e9
        retcode = _solve_retcode_from_error(err)
        if ismissing(retcode)
            retcode = "Exception"
        end
        return (
            ok=false,
            elapsed_s=elapsed_s,
            success=false,
            retcode=String(retcode),
            solver_mode=missing,
            solver_sequence=missing,
            solution=nothing,
            error_text=_perf_error_text(err)
        )
    end
end

function _write_split_rollout_gate_report(
    path::String,
    spec::ProfileSpec,
    gate_df::DataFrame,
    gate_csv_path::String
)
    open(path, "w") do io
        println(io, "# Split-IMEX Rollout Gate (`$(spec.name)` profile)")
        println(io)
        println(io, "- Generated (UTC): $(string(now(UTC)))")
        println(io, "- Cases requested: `$(join(_split_rollout_case_names(), ", "))`")
        println(io, "- Split solvers: `$(join(_split_rollout_solver_variants(), ", "))`")
        println(io, "- Runtime slowdown ceiling: `$(_split_rollout_max_slowdown_ratio())x`")
        println(io, "- Position relative tolerance: `$(_split_rollout_pos_rel_tol())`")
        println(io, "- Velocity relative tolerance: `$(_split_rollout_vel_rel_tol())`")
        println(io, "- Quaternion angle tolerance [rad]: `$(_split_rollout_q_angle_tol_rad())`")
        println(io, "- Angular-rate relative tolerance: `$(_split_rollout_omega_rel_tol())`")
        println(io, "- Trajectory samples per comparison: `$(_split_rollout_sample_count())`")
        println(io, "- Enforce mode: `$(_split_rollout_enforce())`")
        println(io)
        println(io, "- Gate CSV: `$(gate_csv_path)`")
        println(io)
        pass_count = (nrow(gate_df) == 0 || !(:pass_all in names(gate_df))) ? 0 : count(Bool.(gate_df.pass_all))
        println(io, "- Gate pass count: `$(pass_count)/$(nrow(gate_df))`")
        println(io)
        println(io, "| Scenario | Split Solver | Baseline Retcode | Split Retcode | Runtime Ratio | Pos Rel Max | Vel Rel Max | Q Angle Max [rad] | Omega Rel Max | Pass Runtime | Pass Trajectory | Pass All |")
        println(io, "|---|---|---|---|---:|---:|---:|---:|---:|---|---|---|")
        for row in eachrow(gate_df)
            pass_traj = Bool(row.pass_pos) && Bool(row.pass_vel) && Bool(row.pass_q) && Bool(row.pass_omega)
            println(
                io,
                "| $(row.scenario) | $(row.split_solver) | $(row.baseline_retcode) | $(row.split_retcode) | $(_fmt(row.runtime_ratio; digits=4)) | " *
                "$(_fmt(row.pos_rel_max; digits=4)) | $(_fmt(row.vel_rel_max; digits=4)) | $(_fmt(row.q_angle_max_rad; digits=4)) | " *
                "$(_fmt(row.omega_rel_max; digits=4)) | $(row.pass_runtime) | $(pass_traj) | $(row.pass_all) |"
            )
        end
    end
    return nothing
end

function evaluate_split_rollout_gate(
    spec::ProfileSpec,
    cases::Vector{BenchmarkCase},
    outdir::String
)
    requested_names = _split_rollout_case_names()
    split_solvers = _split_rollout_solver_variants()
    case_pool = selected_cases(spec, cases)
    case_by_name = Dict(c.name => c for c in case_pool)
    max_slowdown = _split_rollout_max_slowdown_ratio()
    pos_tol = _split_rollout_pos_rel_tol()
    vel_tol = _split_rollout_vel_rel_tol()
    q_tol = _split_rollout_q_angle_tol_rad()
    omega_tol = _split_rollout_omega_rel_tol()
    sample_count = _split_rollout_sample_count()
    rows = NamedTuple[]

    for scenario_name in requested_names
        if !haskey(case_by_name, scenario_name)
            @warn "[split-rollout] requested scenario '$scenario_name' was not found in profile=$(spec.name); skipping."
            continue
        end
        case = case_by_name[scenario_name]
        baseline_case = _case_with_solver(case; solver_mode_override="auto_stiff", split_imex_solver_override=nothing)

        # Warm up both solver paths so the gate compares runtime behavior, not first-call compilation.
        run_warmup(baseline_case, 1, spec.name)
        for split_solver in split_solvers
            split_case = _case_with_solver(case; solver_mode_override="split_imex", split_imex_solver_override=split_solver)
            run_warmup(split_case, 1, spec.name)

            baseline_run = _run_split_gate_solution(
                baseline_case,
                spec.name;
                solver_mode="auto_stiff",
                split_solver=nothing
            )
            split_run = _run_split_gate_solution(
                split_case,
                spec.name;
                solver_mode="split_imex",
                split_solver=split_solver
            )

            runtime_ratio = (baseline_run.elapsed_s > 0.0) ? split_run.elapsed_s / baseline_run.elapsed_s : Inf
            pass_runtime = baseline_run.success && split_run.success && isfinite(runtime_ratio) && runtime_ratio <= max_slowdown

            pos_rel_max = missing
            vel_rel_max = missing
            q_angle_max_rad = missing
            omega_rel_max = missing
            compared_t_start = missing
            compared_t_end = missing
            compared_samples = missing
            pass_pos = false
            pass_vel = false
            pass_q = !case.args_template.mission_configuration.orientation_sim
            pass_omega = !case.args_template.mission_configuration.orientation_sim

            if baseline_run.success && split_run.success
                metrics = _trajectory_delta_metrics(
                    baseline_run.solution,
                    split_run.solution,
                    length(case.args_template.dynamics_model.spacecraft),
                    case.args_template.mission_configuration.orientation_sim;
                    n_samples=sample_count
                )
                pos_rel_max = metrics.pos_rel_max
                vel_rel_max = metrics.vel_rel_max
                q_angle_max_rad = metrics.q_angle_max_rad
                omega_rel_max = metrics.omega_rel_max
                compared_t_start = metrics.t_start
                compared_t_end = metrics.t_end
                compared_samples = metrics.sample_count
                pass_pos = metrics.pos_rel_max <= pos_tol
                pass_vel = metrics.vel_rel_max <= vel_tol
                if case.args_template.mission_configuration.orientation_sim
                    pass_q = !(metrics.q_angle_max_rad isa Missing) && metrics.q_angle_max_rad <= q_tol
                    pass_omega = !(metrics.omega_rel_max isa Missing) && metrics.omega_rel_max <= omega_tol
                end
            end

            pass_all = pass_runtime && pass_pos && pass_vel && pass_q && pass_omega
            push!(rows, (
                profile=spec.name,
                scenario=scenario_name,
                split_solver=split_solver,
                satellites=length(case.args_template.dynamics_model.spacecraft),
                orientation=case.args_template.mission_configuration.orientation_sim,
                baseline_elapsed_s=baseline_run.elapsed_s,
                split_elapsed_s=split_run.elapsed_s,
                runtime_ratio=runtime_ratio,
                max_slowdown_ratio=max_slowdown,
                pass_runtime=pass_runtime,
                pos_rel_max=pos_rel_max,
                vel_rel_max=vel_rel_max,
                q_angle_max_rad=q_angle_max_rad,
                omega_rel_max=omega_rel_max,
                pos_rel_tol=pos_tol,
                vel_rel_tol=vel_tol,
                q_angle_tol_rad=q_tol,
                omega_rel_tol=omega_tol,
                compared_t_start_s=compared_t_start,
                compared_t_end_s=compared_t_end,
                compared_samples=compared_samples,
                pass_pos=pass_pos,
                pass_vel=pass_vel,
                pass_q=pass_q,
                pass_omega=pass_omega,
                pass_all=pass_all,
                baseline_solver_mode=baseline_run.solver_mode,
                baseline_solver_sequence=baseline_run.solver_sequence,
                baseline_retcode=baseline_run.retcode,
                baseline_error=baseline_run.error_text,
                split_solver_mode=split_run.solver_mode,
                split_solver_sequence=split_run.solver_sequence,
                split_retcode=split_run.retcode,
                split_error=split_run.error_text
            ))
        end
    end

    gate_df = DataFrame(rows)
    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    gate_csv_path = joinpath(outdir, "split_rollout_gate_$(spec.name)_$(stamp).csv")
    gate_report_path = joinpath(outdir, "split_rollout_gate_$(spec.name)_$(stamp).md")
    CSV.write(gate_csv_path, gate_df)
    _write_split_rollout_gate_report(gate_report_path, spec, gate_df, gate_csv_path)

    if _split_rollout_enforce() && nrow(gate_df) > 0 && (:pass_all in names(gate_df)) && any(.!Bool.(gate_df.pass_all))
        failing = gate_df[.!gate_df.pass_all, :]
        summary = join(["$(row.scenario):$(row.split_solver)" for row in eachrow(failing)], ", ")
        error("Split rollout gate failed for $(nrow(failing)) configuration(s): $summary")
    end

    return (df=gate_df, csv_path=gate_csv_path, report_path=gate_report_path)
end

function _write_multirate_rollout_gate_report(
    path::String,
    spec::ProfileSpec,
    gate_df::DataFrame,
    gate_csv_path::String
)
    settings = _multirate_rollout_setting_snapshot()
    open(path, "w") do io
        println(io, "# Multirate Rollout Gate (`$(spec.name)` profile)")
        println(io)
        println(io, "- Generated (UTC): $(string(now(UTC)))")
        println(io, "- Cases requested: `$(join(_multirate_rollout_case_names(), ", "))`")
        println(io, "- Runtime slowdown ceiling: `$(_multirate_rollout_max_slowdown_ratio())x`")
        println(io, "- Position relative tolerance: `$(_multirate_rollout_pos_rel_tol())`")
        println(io, "- Velocity relative tolerance: `$(_multirate_rollout_vel_rel_tol())`")
        println(io, "- Quaternion angle tolerance [rad]: `$(_multirate_rollout_q_angle_tol_rad())`")
        println(io, "- Angular-rate relative tolerance: `$(_multirate_rollout_omega_rel_tol())`")
        println(io, "- Trajectory samples per comparison: `$(_multirate_rollout_sample_count())`")
        println(io, "- Enforce mode: `$(_multirate_rollout_enforce())`")
        println(io, "- Multirate slow solver: `$(settings.slow_solver)`")
        println(io, "- Multirate fast solver: `$(settings.fast_solver)`")
        println(io, "- Multirate slow dt [s]: `$(_fmt(settings.slow_dt_s; digits=4))`")
        println(io, "- Multirate fast substeps: `$(_fmt(settings.fast_substeps; digits=0))`")
        println(io)
        println(io, "- Gate CSV: `$(gate_csv_path)`")
        println(io)
        pass_count = (nrow(gate_df) == 0 || !(:pass_all in names(gate_df))) ? 0 : count(Bool.(gate_df.pass_all))
        println(io, "- Gate pass count: `$(pass_count)/$(nrow(gate_df))`")
        println(io)
        println(io, "| Scenario | Baseline Retcode | Multirate Retcode | Runtime Ratio | Pos Rel Max | Vel Rel Max | Q Angle Max [rad] | Omega Rel Max | Pass Runtime | Pass Trajectory | Pass All |")
        println(io, "|---|---|---|---:|---:|---:|---:|---:|---|---|---|")
        for row in eachrow(gate_df)
            pass_traj = Bool(row.pass_pos) && Bool(row.pass_vel) && Bool(row.pass_q) && Bool(row.pass_omega)
            println(
                io,
                "| $(row.scenario) | $(row.baseline_retcode) | $(row.multirate_retcode) | $(_fmt(row.runtime_ratio; digits=4)) | " *
                "$(_fmt(row.pos_rel_max; digits=4)) | $(_fmt(row.vel_rel_max; digits=4)) | $(_fmt(row.q_angle_max_rad; digits=4)) | " *
                "$(_fmt(row.omega_rel_max; digits=4)) | $(row.pass_runtime) | $(pass_traj) | $(row.pass_all) |"
            )
        end
    end
    return nothing
end

function evaluate_multirate_rollout_gate(
    spec::ProfileSpec,
    cases::Vector{BenchmarkCase},
    outdir::String
)
    requested_names = _multirate_rollout_case_names()
    case_pool = selected_cases(spec, cases)
    case_by_name = Dict(c.name => c for c in case_pool)
    max_slowdown = _multirate_rollout_max_slowdown_ratio()
    pos_tol = _multirate_rollout_pos_rel_tol()
    vel_tol = _multirate_rollout_vel_rel_tol()
    q_tol = _multirate_rollout_q_angle_tol_rad()
    omega_tol = _multirate_rollout_omega_rel_tol()
    sample_count = _multirate_rollout_sample_count()
    settings = _multirate_rollout_setting_snapshot()
    rows = NamedTuple[]

    for scenario_name in requested_names
        if !haskey(case_by_name, scenario_name)
            @warn "[multirate-rollout] requested scenario '$scenario_name' was not found in profile=$(spec.name); skipping."
            continue
        end
        case = case_by_name[scenario_name]
        baseline_case = _case_with_solver(case; solver_mode_override="auto_stiff", split_imex_solver_override=nothing)
        multirate_case = _case_with_solver(case; solver_mode_override="multirate", split_imex_solver_override=nothing)

        # Warm up both solver paths so the gate compares runtime behavior, not first-call compilation.
        run_warmup(baseline_case, 1, spec.name)
        run_warmup(multirate_case, 1, spec.name)

        baseline_run = _run_split_gate_solution(
            baseline_case,
            spec.name;
            solver_mode="auto_stiff",
            split_solver=nothing
        )
        multirate_run = _run_split_gate_solution(
            multirate_case,
            spec.name;
            solver_mode="multirate",
            split_solver=nothing
        )

        runtime_ratio = (baseline_run.elapsed_s > 0.0) ? multirate_run.elapsed_s / baseline_run.elapsed_s : Inf
        pass_runtime = baseline_run.success && multirate_run.success && isfinite(runtime_ratio) && runtime_ratio <= max_slowdown

        pos_rel_max = missing
        vel_rel_max = missing
        q_angle_max_rad = missing
        omega_rel_max = missing
        compared_t_start = missing
        compared_t_end = missing
        compared_samples = missing
        pass_pos = false
        pass_vel = false
        pass_q = !case.args_template.mission_configuration.orientation_sim
        pass_omega = !case.args_template.mission_configuration.orientation_sim

        if baseline_run.success && multirate_run.success
            metrics = _trajectory_delta_metrics(
                baseline_run.solution,
                multirate_run.solution,
                length(case.args_template.dynamics_model.spacecraft),
                case.args_template.mission_configuration.orientation_sim;
                n_samples=sample_count
            )
            pos_rel_max = metrics.pos_rel_max
            vel_rel_max = metrics.vel_rel_max
            q_angle_max_rad = metrics.q_angle_max_rad
            omega_rel_max = metrics.omega_rel_max
            compared_t_start = metrics.t_start
            compared_t_end = metrics.t_end
            compared_samples = metrics.sample_count
            pass_pos = metrics.pos_rel_max <= pos_tol
            pass_vel = metrics.vel_rel_max <= vel_tol
            if case.args_template.mission_configuration.orientation_sim
                pass_q = !(metrics.q_angle_max_rad isa Missing) && metrics.q_angle_max_rad <= q_tol
                pass_omega = !(metrics.omega_rel_max isa Missing) && metrics.omega_rel_max <= omega_tol
            end
        end

        pass_all = pass_runtime && pass_pos && pass_vel && pass_q && pass_omega
        push!(rows, (
            profile=spec.name,
            scenario=scenario_name,
            satellites=length(case.args_template.dynamics_model.spacecraft),
            orientation=case.args_template.mission_configuration.orientation_sim,
            baseline_elapsed_s=baseline_run.elapsed_s,
            multirate_elapsed_s=multirate_run.elapsed_s,
            runtime_ratio=runtime_ratio,
            max_slowdown_ratio=max_slowdown,
            pass_runtime=pass_runtime,
            pos_rel_max=pos_rel_max,
            vel_rel_max=vel_rel_max,
            q_angle_max_rad=q_angle_max_rad,
            omega_rel_max=omega_rel_max,
            pos_rel_tol=pos_tol,
            vel_rel_tol=vel_tol,
            q_angle_tol_rad=q_tol,
            omega_rel_tol=omega_tol,
            compared_t_start_s=compared_t_start,
            compared_t_end_s=compared_t_end,
            compared_samples=compared_samples,
            pass_pos=pass_pos,
            pass_vel=pass_vel,
            pass_q=pass_q,
            pass_omega=pass_omega,
            pass_all=pass_all,
            multirate_slow_solver=settings.slow_solver,
            multirate_fast_solver=settings.fast_solver,
            multirate_slow_dt_s=settings.slow_dt_s,
            multirate_fast_substeps=settings.fast_substeps,
            baseline_solver_mode=baseline_run.solver_mode,
            baseline_solver_sequence=baseline_run.solver_sequence,
            baseline_retcode=baseline_run.retcode,
            baseline_error=baseline_run.error_text,
            multirate_solver_mode=multirate_run.solver_mode,
            multirate_solver_sequence=multirate_run.solver_sequence,
            multirate_retcode=multirate_run.retcode,
            multirate_error=multirate_run.error_text
        ))
    end

    gate_df = DataFrame(rows)
    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    gate_csv_path = joinpath(outdir, "multirate_rollout_gate_$(spec.name)_$(stamp).csv")
    gate_report_path = joinpath(outdir, "multirate_rollout_gate_$(spec.name)_$(stamp).md")
    CSV.write(gate_csv_path, gate_df)
    _write_multirate_rollout_gate_report(gate_report_path, spec, gate_df, gate_csv_path)

    if _multirate_rollout_enforce() && nrow(gate_df) > 0 && (:pass_all in names(gate_df)) && any(.!Bool.(gate_df.pass_all))
        failing = gate_df[.!gate_df.pass_all, :]
        summary = join([String(row.scenario) for row in eachrow(failing)], ", ")
        error("Multirate rollout gate failed for $(nrow(failing)) configuration(s): $summary")
    end

    return (df=gate_df, csv_path=gate_csv_path, report_path=gate_report_path)
end

@inline _plot_label(name::AbstractString) = replace(name, "_" => " ")
@inline _plot_axis_label(name::AbstractString) = replace(name, "_" => "\n")
@inline _plot_number(v) = v isa Missing ? NaN : Float64(v)

function _plot_wrapped_label(name::AbstractString; width::Int=16, max_lines::Int=3)::String
    words = split(_plot_label(name))
    isempty(words) && return ""
    lines = String[]
    current = ""
    for word in words
        if isempty(current)
            current = word
        elseif ncodeunits(current) + 1 + ncodeunits(word) <= width
            current *= " " * word
        else
            push!(lines, current)
            current = word
        end
    end
    push!(lines, current)
    if length(lines) > max_lines
        head = lines[1:(max_lines - 1)]
        tail = join(lines[max_lines:end], " ")
        push!(head, tail)
        lines = head
    end
    return join(lines, "\n")
end

@inline function _plot_ready()::Bool
    return myid() == 1 && isdefined(Main, :Plots)
end

const _runtime_plot_theme_applied = Ref(false)

function _ensure_runtime_plot_theme!()
    !_plot_ready() && return nothing
    _runtime_plot_theme_applied[] && return nothing

    Plots.theme(:ggplot2)
    Plots.default(
        dpi=220,
        lw=2,
        ms=5,
        markerstrokewidth=1.3,
        markerstrokecolor=:black,
        titlefont=Plots.font(24),
        guidefont=Plots.font(16),
        tickfont=Plots.font(12),
        legend_font=Plots.font(11),
        legend_background_color=:white,
        legend_foreground_color=:black,
        gridalpha=0.24,
        minorgrid=false,
        framestyle=:box
    )
    _runtime_plot_theme_applied[] = true
    return nothing
end

@inline function _plot_margins(;
    size::Tuple{Int, Int}=(2200, 1100),
    left_mm::Int=18,
    right_mm::Int=18,
    top_mm::Int=8,
    bottom_mm::Int=20,
    legend=false,
    xrotation::Real=0
)
    return (
        size=size,
        left_margin=left_mm * Plots.mm,
        right_margin=right_mm * Plots.mm,
        top_margin=top_mm * Plots.mm,
        bottom_margin=bottom_mm * Plots.mm,
        legend=legend,
        xrotation=xrotation,
        framestyle=:box,
        gridalpha=0.24,
        legend_background_color=:white,
        legend_foreground_color=:black
    )
end

function _plot_metric_pairs(df::DataFrame, label_col::Symbol, metric_col::Symbol; axis_labels::Bool=false)
    labels = String[]
    values = Float64[]
    for row in eachrow(df)
        value = row[metric_col]
        if !(value isa Missing)
            f = Float64(value)
            if isfinite(f)
                raw = string(row[label_col])
                push!(labels, axis_labels ? _plot_axis_label(raw) : _plot_label(raw))
                push!(values, f)
            end
        end
    end
    return labels, values
end

@inline function _has_row_fields(row, fields::Vector{Symbol})::Bool
    for field in fields
        if row[field] isa Missing
            return false
        end
    end
    return true
end

function _save_runtime_plot!(
    artifacts::Vector{String},
    plt,
    outdir::String,
    basename::String,
    spec::ProfileSpec,
    stamp::String
)
    path = joinpath(outdir, "$(basename)_$(spec.name)_$(stamp).png")
    try
        Plots.savefig(plt, path)
        push!(artifacts, path)
    catch err
        @warn "[perf] failed to save plot $basename: $(_perf_error_text(err))"
    end
    return nothing
end

@inline function _paper_figure_pack_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_PAPER_FIGURE_PACK", true)
end

@inline function _paper_ladder_outdir()::String
    raw = strip(get(ENV, "SPACEAGORA_PERF_PAPER_LADDER_OUTDIR", joinpath(DEFAULT_OUTPUT_DIR, "smart_parallel_ladder")))
    return normpath(isabspath(raw) ? raw : joinpath(REPO_ROOT, raw))
end

@inline function _paper_cross_machine_outdir()::String
    raw = strip(get(ENV, "SPACEAGORA_PERF_PAPER_CROSS_MACHINE_OUTDIR", joinpath(DEFAULT_OUTPUT_DIR, "smart_parallel_ladder_cross_machine")))
    return normpath(isabspath(raw) ? raw : joinpath(REPO_ROOT, raw))
end

function _latest_profile_artifact_optional(
    outdir::String,
    prefix::String,
    profile::String;
    suffix::String=".csv"
)::Union{Nothing, String}
    isdir(outdir) || return nothing
    pattern = "$(prefix)_$(profile)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

function _read_optional_dataframe(path::Union{Nothing, String}, label::String)::Union{Nothing, DataFrame}
    path === nothing && return nothing
    try
        return CSV.read(path, DataFrame)
    catch err
        @warn "[perf] Failed to read optional paper-figure artifact; skipping." label=label path=path exception=(err, catch_backtrace())
        return nothing
    end
end

function _paper_figure_external_data(spec::ProfileSpec)
    ladder_outdir = _paper_ladder_outdir()
    cross_outdir = _paper_cross_machine_outdir()

    layer_attr_speedup_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_layer_attribution_speedup",
        spec.name
    )
    deep_accuracy_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_deep_accuracy_parity",
        spec.name
    )
    montecarlo_parity_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_montecarlo_distribution_parity",
        spec.name
    )
    route_mix_path = _latest_profile_artifact_optional(
        ladder_outdir,
        "smart_parallel_ladder_route_mix",
        spec.name
    )
    cross_speedup_summary_path = _latest_profile_artifact_optional(
        cross_outdir,
        "smart_parallel_ladder_cross_machine_speedup_summary",
        spec.name
    )
    cross_adaptive_regret_summary_path = _latest_profile_artifact_optional(
        cross_outdir,
        "smart_parallel_ladder_cross_machine_adaptive_regret_summary",
        spec.name
    )
    cross_route_mix_summary_path = _latest_profile_artifact_optional(
        cross_outdir,
        "smart_parallel_ladder_cross_machine_route_mix_summary",
        spec.name
    )

    return (
        layer_attribution_speedup_df=_read_optional_dataframe(layer_attr_speedup_path, "layer_attribution_speedup"),
        deep_accuracy_df=_read_optional_dataframe(deep_accuracy_path, "deep_accuracy_parity"),
        montecarlo_parity_df=_read_optional_dataframe(montecarlo_parity_path, "montecarlo_distribution_parity"),
        route_mix_df=_read_optional_dataframe(route_mix_path, "route_mix"),
        cross_speedup_summary_df=_read_optional_dataframe(cross_speedup_summary_path, "cross_machine_speedup_summary"),
        cross_adaptive_regret_summary_df=_read_optional_dataframe(cross_adaptive_regret_summary_path, "cross_machine_adaptive_regret_summary"),
        cross_route_mix_summary_df=_read_optional_dataframe(cross_route_mix_summary_path, "cross_machine_route_mix_summary"),
        layer_attribution_speedup_path=layer_attr_speedup_path,
        deep_accuracy_path=deep_accuracy_path,
        montecarlo_parity_path=montecarlo_parity_path,
        route_mix_path=route_mix_path,
        cross_speedup_summary_path=cross_speedup_summary_path,
        cross_adaptive_regret_summary_path=cross_adaptive_regret_summary_path,
        cross_route_mix_summary_path=cross_route_mix_summary_path
    )
end

function _sorted_orbit_groups(df::DataFrame, metric::Symbol)
    multiplier_col = :mission_time_multiplier in names(df) ? :mission_time_multiplier : :orbit_count
    groups = collect(groupby(df, :scenario))
    sort!(groups; by=g -> begin
        local_df = DataFrame(g)
        sort!(local_df, multiplier_col)
        value = local_df[1, metric]
        return value isa Missing ? Inf : -Float64(value)
    end)
    return groups
end

function generate_runtime_plots(
    outdir::String,
    spec::ProfileSpec,
    stamp::String,
    raw_df::DataFrame,
    summary_df::DataFrame,
    orbit_summary_df::DataFrame,
    entry_duration_summary_df::DataFrame=DataFrame();
    split_gate_df::Union{Nothing, DataFrame}=nothing,
    multirate_gate_df::Union{Nothing, DataFrame}=nothing,
    inner_hint_layer_df::Union{Nothing, DataFrame}=nothing,
    density_backend_breakdown_df::Union{Nothing, DataFrame}=nothing
)::Vector{String}
    !_plot_ready() && return String[]
    _ensure_runtime_plot_theme!()

    plot_artifacts = String[]
    success_summary = summary_df[summary_df.samples_success .> 0, :]
    paper_external = _paper_figure_pack_enabled() ? _paper_figure_external_data(spec) : nothing

    # 1) Mean and p90 runtime per scenario.
    totals_df = success_summary[.!ismissing.(success_summary.total_time_mean_s), :]
    if nrow(totals_df) > 0
        labels = _plot_axis_label.(String.(totals_df.scenario))
        means = Float64.(totals_df.total_time_mean_s)
        p90_vals = [_plot_number(v) for v in totals_df.total_time_p90_s]
        valid_p90 = findall(isfinite, p90_vals)

        plt = Plots.bar(
            labels,
            means;
            label="Mean total time",
            color="#5b8fb9",
            title="Scenario Runtime: Mean + P90 (bars sorted by mean)",
            xlabel="Scenario",
            ylabel="Wall Time [s]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        if !isempty(valid_p90)
            Plots.scatter!(
                plt,
                labels[valid_p90],
                p90_vals[valid_p90];
                marker=:diamond,
                color=:black,
                label="P90 total time"
            )
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_totals", spec, stamp)
    end

    # 2) Relative speedup vs selected baseline scenario.
    speedup_df = success_summary[.!ismissing.(success_summary.speedup_vs_baseline), :]
    if nrow(speedup_df) > 0
        labels = _plot_axis_label.(String.(speedup_df.scenario))
        speedups = Float64.(speedup_df.speedup_vs_baseline)
        plt = Plots.bar(
            labels,
            speedups;
            label="Speedup",
            color="#4f9d69",
            title="Relative Speedup (higher is faster)",
            xlabel="Scenario",
            ylabel="Speedup vs $(PERF_BASELINE_SCENARIO) [x]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        Plots.hline!(plt, [1.0]; color=:black, linestyle=:dash, label="Baseline = 1x")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_speedup", spec, stamp)
    end

    # 3) Robustness-adjusted performance (all attempts) + failure rate.
    if (:penalized_expected_wall_time_s in Symbol.(names(summary_df))) &&
       (:failure_rate in Symbol.(names(summary_df)))
        robust_df = summary_df[.!ismissing.(summary_df.penalized_expected_wall_time_s), :]
        if nrow(robust_df) > 0
            sort!(robust_df, :penalized_expected_wall_time_s)
            labels = _plot_axis_label.(String.(robust_df.scenario))
            penalized = Float64.(robust_df.penalized_expected_wall_time_s)
            failure_pct = [
                v isa Missing ? 0.0 : (100.0 * Float64(v))
                for v in robust_df.failure_rate
            ]
            p1 = Plots.bar(
                labels,
                penalized;
                color="#3f7fb3",
                label=false,
                title="Robustness-Adjusted Performance (All Attempts)",
                xlabel="Scenario",
                ylabel="Penalized Expected Wall Time [s / requested run]",
                _plot_margins(size=(2500, 760), bottom_mm=92, right_mm=32)...
            )
            p2 = Plots.bar(
                labels,
                failure_pct;
                color="#bc4749",
                label=false,
                title="Failure Rate by Scenario (All Attempts)",
                xlabel="Scenario",
                ylabel="Failure Rate [%]",
                _plot_margins(size=(2500, 760), bottom_mm=92, right_mm=32)...
            )
            plt = Plots.plot(
                p1,
                p2;
                layout=(2, 1),
                size=(2500, 1600),
                left_margin=20 * Plots.mm,
                right_margin=18 * Plots.mm
            )
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_robustness_adjusted", spec, stamp)
        end
    end

    # 4) Runtime variability by scenario.
    variability_df = success_summary[
        .!ismissing.(success_summary.total_time_min_s) .&
        .!ismissing.(success_summary.total_time_mean_s) .&
        .!ismissing.(success_summary.total_time_p90_s) .&
        .!ismissing.(success_summary.total_time_max_s), :
    ]
    if nrow(variability_df) > 0
        labels = _plot_axis_label.(String.(variability_df.scenario))
        x = collect(1:length(labels))
        mins = Float64.(variability_df.total_time_min_s)
        means = Float64.(variability_df.total_time_mean_s)
        p90s = Float64.(variability_df.total_time_p90_s)
        maxs = Float64.(variability_df.total_time_max_s)
        plt = Plots.plot(
            x,
            means;
            label="Mean",
            color="#3d83b8",
            marker=:circle,
            xticks=(x, labels),
            title="Runtime Variability by Scenario (Min/Mean/P90/Max)",
            xlabel="Scenario",
            ylabel="Wall Time [s]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        for i in eachindex(x)
            Plots.plot!(
                plt,
                [x[i], x[i]],
                [mins[i], maxs[i]];
                color=:gray40,
                linewidth=2,
                label=i == 1 ? "Min-Max range" : ""
            )
        end
        Plots.scatter!(plt, x, p90s; marker=:diamond, color=:black, label="P90")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_variability", spec, stamp)
    end

    # 5) Configuration copy + solve breakdown.
    labels_breakdown = String[]
    copy_vals = Float64[]
    solve_vals = Float64[]
    for row in eachrow(success_summary)
        if _has_row_fields(row, [:copy_time_mean_s, :solve_time_mean_s])
            push!(labels_breakdown, _plot_axis_label(String(row.scenario)))
            push!(copy_vals, Float64(row.copy_time_mean_s))
            push!(solve_vals, Float64(row.solve_time_mean_s))
        end
    end
    if !isempty(labels_breakdown)
        plt = Plots.bar(
            labels_breakdown,
            hcat(copy_vals, solve_vals);
            label=["Copy time" "Solve time"],
            color=["#9fb3c8" "#2a9d8f"],
            title="Runtime Breakdown (Configuration Copy + Solve)",
            xlabel="Scenario",
            ylabel="Wall Time [s]",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_breakdown_copy_solve", spec, stamp)
    end

    # 6) Allocation footprint by scenario (memory + call count).
    alloc_df = success_summary[
        .!ismissing.(success_summary.total_bytes_mean_mb) .&
        .!ismissing.(success_summary.solve_alloc_mean), :
    ]
    if nrow(alloc_df) > 0
        labels = _plot_axis_label.(String.(alloc_df.scenario))
        x = collect(1:length(labels))
        bytes_mb = Float64.(alloc_df.total_bytes_mean_mb)
        alloc_million = Float64.(alloc_df.solve_alloc_mean) ./ 1e6
        p1 = Plots.bar(
            x,
            bytes_mb;
            color="#d17a4f",
            label=false,
            xticks=(x, fill("", length(x))),
            ylabel="Allocated Memory [MB]",
            title="Allocation Footprint by Scenario",
            _plot_margins(size=(2500, 780), bottom_mm=10)...
        )
        p2 = Plots.bar(
            x,
            alloc_million;
            color="#6999a1",
            label=false,
            xticks=(x, labels),
            xlabel="Scenario",
            ylabel="Allocation Calls [million]",
            _plot_margins(size=(2500, 980), bottom_mm=88)...
        )
        plt = Plots.plot(
            p1,
            p2;
            layout=(2, 1),
            size=(2500, 1760),
            left_margin=20 * Plots.mm,
            right_margin=20 * Plots.mm
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_memory_alloc", spec, stamp)
    end

    # 7) Integrator workload and rejection pressure.
    solver_df = success_summary[
        .!ismissing.(success_summary.accepted_steps_mean) .&
        .!ismissing.(success_summary.saved_points_mean), :
    ]
    if nrow(solver_df) > 0
        labels = _plot_axis_label.(String.(solver_df.scenario))
        x = collect(1:length(labels))
        accepted = Float64.(solver_df.accepted_steps_mean)
        saved = Float64.(solver_df.saved_points_mean)
        rejected = [v isa Missing ? 0.0 : Float64(v) for v in solver_df.rejected_steps_mean]
        rejection_ratio = [acc + rej <= 0.0 ? 0.0 : 100.0 * rej / (acc + rej) for (acc, rej) in zip(accepted, rejected)]
        p1 = Plots.plot(
            x,
            accepted;
            label="Accepted steps",
            color="#2e7d32",
            marker=:circle,
            xticks=(x, fill("", length(x))),
            ylabel="Steps / Saved Points",
            title="Integrator Workload per Scenario",
            _plot_margins(size=(2500, 780), bottom_mm=10, legend=:outertopright)...
        )
        Plots.plot!(p1, x, saved; color="#375fd2", marker=:diamond, label="Saved points")
        p2 = Plots.bar(
            x,
            rejection_ratio;
            color="#cc6666",
            label=false,
            xticks=(x, labels),
            title="Solver Rejection Pressure",
            xlabel="Scenario",
            ylabel="Rejected Step Ratio [%]",
            _plot_margins(size=(2500, 980), bottom_mm=88)...
        )
        plt = Plots.plot(
            p1,
            p2;
            layout=(2, 1),
            size=(2500, 1760),
            left_margin=20 * Plots.mm,
            right_margin=20 * Plots.mm
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_solver_workload", spec, stamp)
    end

    # 8) Throughput ranking with per-satellite markers.
    throughput_df = success_summary[.!ismissing.(success_summary.sim_seconds_per_wall_second_mean), :]
    if nrow(throughput_df) > 0
        sort!(throughput_df, :sim_seconds_per_wall_second_mean, rev=true)
        labels = _plot_axis_label.(String.(throughput_df.scenario))
        x = collect(1:length(labels))
        global_tp = Float64.(throughput_df.sim_seconds_per_wall_second_mean)
        sat_tp = [v isa Missing ? NaN : Float64(v) for v in throughput_df.satellite_sim_seconds_per_wall_second_mean]
        plt = Plots.bar(
            x,
            global_tp;
            color="#4f9d69",
            label="Global throughput",
            xticks=(x, labels),
            title="Simulation Throughput Ranking",
            xlabel="Scenario",
            ylabel="Sim Seconds / Wall Second",
            _plot_margins(size=(2500, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        valid_sat = findall(isfinite, sat_tp)
        if !isempty(valid_sat)
            Plots.scatter!(
                plt,
                x[valid_sat],
                sat_tp[valid_sat];
                marker=:utriangle,
                color=:black,
                label="Per-satellite throughput"
            )
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_throughput", spec, stamp)
    end

    # 9) Satellite-count scaling (measured vs ideal linear).
    sat_df = success_summary[
        (success_summary.category .== "satellite_scaling") .&
        .!ismissing.(success_summary.satellites) .&
        .!ismissing.(success_summary.total_time_mean_s), :
    ]
    if nrow(sat_df) > 0
        sort!(sat_df, :satellites)
        sat_count = Int.(sat_df.satellites)
        measured = Float64.(sat_df.total_time_mean_s)
        base_runtime = measured[1]
        ideal = base_runtime .* (sat_count ./ sat_count[1])
        plt = Plots.plot(
            sat_count,
            measured;
            marker=:circle,
            color="#3f7fb3",
            label="Measured runtime",
            title="Satellite-Count Scaling (Measured vs Ideal Linear)",
            xlabel="Number of Satellites",
            ylabel="Mean Runtime [s]",
            _plot_margins(size=(2200, 1200), bottom_mm=20, right_mm=35, legend=:topleft)...
        )
        Plots.plot!(
            plt,
            sat_count,
            ideal;
            marker=:diamond,
            color=:black,
            linestyle=:dash,
            label="Ideal linear from smallest satellite-count case"
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_satellite_scaling", spec, stamp)
    end

    # 10) Dynamics fidelity ladder (absolute + relative).
    fidelity_order = [
        "single_j2",
        "single_nbody_sun_moon",
        "single_harmonics_l20",
        "single_harmonics_l50",
    ]
    fidelity_labels = String[]
    fidelity_values = Float64[]
    for scenario in fidelity_order
        idx = findfirst(==(scenario), success_summary.scenario)
        if idx !== nothing
            value = success_summary.total_time_mean_s[idx]
            if !(value isa Missing)
                if scenario == "single_j2"
                    push!(fidelity_labels, "J2")
                elseif scenario == "single_nbody_sun_moon"
                    push!(fidelity_labels, "NBody\nSun+Moon")
                elseif scenario == "single_harmonics_l20"
                    push!(fidelity_labels, "Harmonics\nL20")
                elseif scenario == "single_harmonics_l50"
                    push!(fidelity_labels, "Harmonics\nL50")
                else
                    push!(fidelity_labels, _plot_label(scenario))
                end
                push!(fidelity_values, Float64(value))
            end
        end
    end
    if length(fidelity_values) >= 2
        x = collect(1:length(fidelity_values))
        baseline = fidelity_values[1]
        relative = baseline <= 0.0 ? fill(NaN, length(fidelity_values)) : fidelity_values ./ baseline
        p_abs = Plots.bar(
            x,
            fidelity_values;
            color="#7668c7",
            label=false,
            xticks=(x, fill("", length(x))),
            ylabel="Mean Runtime [s]",
            title="Dynamics Fidelity Ladder (Absolute Runtime)",
            _plot_margins(size=(2200, 760), bottom_mm=10)...
        )
        p_rel = Plots.bar(
            x,
            relative;
            color="#d67c1c",
            label=false,
            xticks=(x, fidelity_labels),
            xlabel="Dynamics Model",
            ylabel="Runtime / Baseline [x]",
            title="Dynamics Fidelity Ladder (Relative Runtime)",
            _plot_margins(size=(2200, 940), bottom_mm=56)...
        )
        Plots.hline!(p_rel, [1.0]; color=:black, linestyle=:dash, label=false)
        plt = Plots.plot(
            p_abs,
            p_rel;
            layout=(2, 1),
            size=(2200, 1700),
            left_margin=20 * Plots.mm,
            right_margin=18 * Plots.mm
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_fidelity_ladder", spec, stamp)
    end

    # 11) Monte Carlo runtime distribution with mean and p90 lines.
    mc_df = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
    mc_times = [Float64(v) for v in mc_df.total_time_s if !(v isa Missing)]
    if !isempty(mc_times)
        mc_mean = mean(mc_times)
        mc_p90 = quantile(mc_times, 0.9)
        plt = Plots.histogram(
            mc_times;
            bins=min(20, max(5, round(Int, sqrt(length(mc_times))))),
            color="#7b63c6",
            alpha=0.65,
            label="Samples",
            title="Monte Carlo Runtime Distribution",
            xlabel="Total Wall Time [s]",
            ylabel="Count",
            _plot_margins(size=(2200, 1300), right_mm=52, legend=:outertopright)...
        )
        Plots.vline!(plt, [mc_mean]; color=:red, linewidth=2, label="Mean")
        Plots.vline!(plt, [mc_p90]; color=:black, linewidth=2, linestyle=:dash, label="P90")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_montecarlo_hist", spec, stamp)
    end

    # 12) Monte Carlo seed trace with mean and p90 lines.
    mc_seed_df = mc_df[.!ismissing.(mc_df.seed), :]
    if nrow(mc_seed_df) > 0
        sort!(mc_seed_df, :seed)
        seeds = Int.(mc_seed_df.seed)
        totals = Float64.(mc_seed_df.total_time_s)
        seed_mean = mean(totals)
        seed_p90 = quantile(totals, 0.9)
        plt = Plots.plot(
            seeds,
            totals;
            marker=:circle,
            color="#356a97",
            label="Seed runtime",
            title="Monte Carlo Runtime by Seed",
            xlabel="Monte Carlo Seed",
            ylabel="Total Wall Time [s]",
            _plot_margins(size=(2200, 1300), right_mm=52, legend=:outertopright)...
        )
        Plots.hline!(plt, [seed_mean]; color=:red, linewidth=2, label="Mean")
        Plots.hline!(plt, [seed_p90]; color=:black, linewidth=2, linestyle=:dash, label="P90")
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_montecarlo_seed_trace", spec, stamp)
    end

    sweep_multiplier_col = :mission_time_multiplier in names(orbit_summary_df) ? :mission_time_multiplier : :orbit_count
    orbit_valid = orbit_summary_df[
        (orbit_summary_df.samples_success .> 0) .&
        .!ismissing.(orbit_summary_df[!, sweep_multiplier_col]), :
    ]

    # 13) Mission-time sweep runtime scaling.
    orbit_scaling_df = orbit_valid[.!ismissing.(orbit_valid.total_time_mean_s), :]
    if nrow(orbit_scaling_df) > 0
        multiplier_col = :mission_time_multiplier in names(orbit_scaling_df) ? :mission_time_multiplier : :orbit_count
        time_per_unit_col = :time_per_baseline_period_mean_s in names(orbit_scaling_df) ? :time_per_baseline_period_mean_s : :time_per_orbit_mean_s
        plt = Plots.plot(;
            title="Mission-Time Sweep Runtime Scaling",
            xlabel="Mission-Time Multiplier [x baseline period]",
            ylabel="Mean Time per Baseline-Period Unit [s]",
            _plot_margins(size=(2300, 1200), right_mm=72, legend=:outerright)...
        )
        for grp in _sorted_orbit_groups(orbit_scaling_df, :total_time_mean_s)
            local_df = DataFrame(grp)
            sort!(local_df, multiplier_col)
            x = Int.(local_df[!, multiplier_col])
            y = Float64.(local_df[!, time_per_unit_col])
            Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_scaling", spec, stamp)
    end

    # 14) Mission-time sweep efficiency scaling.
    orbit_eff_df = orbit_valid[.!ismissing.(orbit_valid.orbits_per_wall_second_mean), :]
    if nrow(orbit_eff_df) > 0
        multiplier_col = :mission_time_multiplier in names(orbit_eff_df) ? :mission_time_multiplier : :orbit_count
        throughput_col = :baseline_periods_per_wall_second_mean in names(orbit_eff_df) ? :baseline_periods_per_wall_second_mean : :orbits_per_wall_second_mean
        plt = Plots.plot(;
            title="Mission-Time Sweep Efficiency Scaling",
            xlabel="Mission-Time Multiplier [x baseline period]",
            ylabel="Baseline-Period Units / Wall-sec",
            _plot_margins(size=(2300, 1200), right_mm=72, legend=:outerright)...
        )
        for grp in _sorted_orbit_groups(orbit_eff_df, throughput_col)
            local_df = DataFrame(grp)
            sort!(local_df, multiplier_col)
            x = Int.(local_df[!, multiplier_col])
            y = Float64.(local_df[!, throughput_col])
            Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
        end
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_efficiency", spec, stamp)
    end

    # 15) Mission-time sweep time heatmap.
    heat_df = orbit_valid[.!ismissing.(orbit_valid.time_per_orbit_mean_s), :]
    if nrow(heat_df) > 0
        multiplier_col = :mission_time_multiplier in names(heat_df) ? :mission_time_multiplier : :orbit_count
        heat_value_col = :time_per_baseline_period_mean_s in names(heat_df) ? :time_per_baseline_period_mean_s : :time_per_orbit_mean_s
        scenario_names = unique(String.(heat_df.scenario))
        scenario_order = sort(scenario_names; by=sc -> begin
            vals = [
                Float64(v)
                for v in heat_df[heat_df.scenario .== sc, heat_value_col]
                if !(v isa Missing)
            ]
            return isempty(vals) ? Inf : -mean(vals)
        end)
        multipliers = sort(unique(Int.(heat_df[!, multiplier_col])))
        z = fill(NaN, length(scenario_order), length(multipliers))
        for row in eachrow(heat_df)
            si = findfirst(==(String(row.scenario)), scenario_order)
            oi = findfirst(==(Int(row[ multiplier_col ])), multipliers)
            if si !== nothing && oi !== nothing
                z[si, oi] = Float64(row[heat_value_col])
            end
        end
        plt = Plots.heatmap(
            multipliers,
            1:length(scenario_order),
            z;
            color=Plots.cgrad([:lightsteelblue1, :mediumpurple3, "#3a0a2a"]),
            colorbar_title="s / baseline-period unit",
            colorbar=true,
            yticks=(1:length(scenario_order), _plot_wrapped_label.(scenario_order)),
            xlabel="Mission-Time Multiplier [x baseline period]",
            ylabel="Scenario",
            title="Mission-Time Sweep Heatmap [s / baseline-period unit]",
            _plot_margins(size=(2300, 1400), left_mm=72, right_mm=52, bottom_mm=24)...
        )
        _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_per_orbit_heatmap", spec, stamp)
    end

    # 16) Entry-duration sweep trends (separate from per-orbit mission-time sweep).
    if nrow(entry_duration_summary_df) > 0 &&
       (:entry_run_role in names(entry_duration_summary_df)) &&
       (:entry_atmospheric_interface_count in names(entry_duration_summary_df))
        entry_measured = entry_duration_summary_df[
            (entry_duration_summary_df.entry_run_role .== "measured") .&
            .!ismissing.(entry_duration_summary_df.entry_atmospheric_interface_count), :
        ]
        if nrow(entry_measured) > 0
            if :passage_duration_mean_s in names(entry_measured)
                passage_df = entry_measured[.!ismissing.(entry_measured.passage_duration_mean_s), :]
                if nrow(passage_df) > 0
                    plt = Plots.plot(;
                        title="Entry-Duration Sweep: Passage Duration",
                        xlabel="Atmospheric-Interface Count",
                        ylabel="Passage Duration [s]",
                        _plot_margins(size=(2200, 1200), right_mm=72, legend=:outerright)...
                    )
                    for grp in groupby(passage_df, :scenario)
                        local_df = DataFrame(grp)
                        sort!(local_df, :entry_atmospheric_interface_count)
                        x = Int.(local_df.entry_atmospheric_interface_count)
                        y = Float64.(local_df.passage_duration_mean_s)
                        Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
                    end
                    _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_entry_duration_passage", spec, stamp)
                end
            end

            if :wall_time_per_passage_mean_s in names(entry_measured)
                wall_df = entry_measured[.!ismissing.(entry_measured.wall_time_per_passage_mean_s), :]
                if nrow(wall_df) > 0
                    plt = Plots.plot(;
                        title="Entry-Duration Sweep: Wall Time per Entry Passage",
                        xlabel="Atmospheric-Interface Count",
                        ylabel="Wall Time per Passage [s]",
                        _plot_margins(size=(2200, 1200), right_mm=72, legend=:outerright)...
                    )
                    for grp in groupby(wall_df, :scenario)
                        local_df = DataFrame(grp)
                        sort!(local_df, :entry_atmospheric_interface_count)
                        x = Int.(local_df.entry_atmospheric_interface_count)
                        y = Float64.(local_df.wall_time_per_passage_mean_s)
                        Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
                    end
                    _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_entry_duration_wall_per_passage", spec, stamp)
                end
            end

            if :event_time_abs_error_mean_s in names(entry_measured)
                error_df = entry_measured[.!ismissing.(entry_measured.event_time_abs_error_mean_s), :]
                if nrow(error_df) > 0
                    plt = Plots.plot(;
                        title="Entry-Duration Sweep: Event-Time Absolute Error vs Serial Reference",
                        xlabel="Atmospheric-Interface Count",
                        ylabel="|Event-Time Error| [s]",
                        _plot_margins(size=(2200, 1200), right_mm=72, legend=:outerright)...
                    )
                    for grp in groupby(error_df, :scenario)
                        local_df = DataFrame(grp)
                        sort!(local_df, :entry_atmospheric_interface_count)
                        x = Int.(local_df.entry_atmospheric_interface_count)
                        y = Float64.(local_df.event_time_abs_error_mean_s)
                        Plots.plot!(plt, x, y; marker=:circle, label=String(local_df.scenario[1]))
                    end
                    _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_entry_duration_event_time_error", spec, stamp)
                end
            end
        end
    end

    # 16) SPICE call budget by scenario (N-body + SRP + planet frame).
    spice_df = success_summary[
        .!ismissing.(success_summary.spice_calls_total_mean) .&
        .!ismissing.(success_summary.total_time_mean_s), :
    ]
    if nrow(spice_df) > 0
        labels = _plot_axis_label.(String.(spice_df.scenario))
        nbody_calls = [v isa Missing ? 0.0 : Float64(v) for v in spice_df.nbody_spkpos_total_calls_mean]
        srp_calls = [v isa Missing ? 0.0 : Float64(v) for v in spice_df.srp_spkpos_total_calls_mean]
        pxform_calls = [v isa Missing ? 0.0 : Float64(v) for v in spice_df.planet_pxform_total_calls_mean]
        total_calls = nbody_calls .+ srp_calls .+ pxform_calls
        calls_per_wall_s = [
            (tt <= 0.0) ? NaN : tc / tt
            for (tc, tt) in zip(total_calls, Float64.(spice_df.total_time_mean_s))
        ]
        if any(>(0.0), total_calls)
            p1 = Plots.bar(
                labels,
                hcat(nbody_calls, srp_calls, pxform_calls);
                label=["N-body spkpos" "SRP spkpos" "Planet pxform"],
                color=["#2f80ed" "#56b870" "#d17a4f"],
                bar_position=:stack,
                ylabel="Mean SPICE Calls / Successful Run",
                title="SPICE Call Budget by Scenario",
                _plot_margins(size=(2500, 760), bottom_mm=10, right_mm=62, legend=:outertopright)...
            )
            p2 = Plots.bar(
                labels,
                calls_per_wall_s;
                color="#7b63c6",
                label=false,
                xlabel="Scenario",
                ylabel="Calls / Wall-sec",
                title="SPICE Call Rate by Scenario",
                _plot_margins(size=(2500, 940), bottom_mm=88, right_mm=42)...
            )
            plt = Plots.plot(
                p1,
                p2;
                layout=(2, 1),
                size=(2500, 1700),
                left_margin=20 * Plots.mm,
                right_margin=20 * Plots.mm
            )
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_spice_budget", spec, stamp)
        end
    end

    # 17) Paper pack: adaptivity behavior (route mix + confidence/regret + best-fixed regret).
    if _paper_figure_pack_enabled()
        panels = Any[]

        route_mix_df = (paper_external === nothing) ? nothing : paper_external.route_mix_df
        if !(route_mix_df === nothing) &&
           nrow(route_mix_df) > 0 &&
           (:none_pct in names(route_mix_df)) &&
           (:threads_pct in names(route_mix_df)) &&
           (:process_pct in names(route_mix_df))
            local_df = copy(route_mix_df)
            if :mode in names(local_df)
                mask = [_is_adaptive_mode_token(v) for v in local_df.mode]
                if any(mask)
                    local_df = local_df[mask, :]
                end
            end
            if nrow(local_df) > 0
                labels = if :rung in names(local_df)
                    _plot_axis_label.(String.(local_df.rung))
                elseif :mode in names(local_df)
                    _plot_axis_label.(String.(local_df.mode))
                else
                    _plot_axis_label.(string.("adaptive_", collect(1:nrow(local_df))))
                end
                none_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.none_pct]
                threads_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.threads_pct]
                process_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.process_pct]
                p_route = Plots.bar(
                    labels,
                    hcat(none_pct, threads_pct, process_pct);
                    label=["none" "threads" "process"],
                    bar_position=:stack,
                    color=["#9aa4b2" "#3f7fb3" "#2f8f5b"],
                    title="Adaptive Route-Choice Distribution",
                    xlabel="Adaptive Mode / Rung",
                    ylabel="Route Share [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_route)
            end
        elseif nrow(raw_df) > 0 && (:outer_route in names(raw_df))
            counts = Dict("none" => 0, "threads" => 0, "process" => 0, "other" => 0)
            for v in raw_df.outer_route
                token = lowercase(strip(String(v)))
                if haskey(counts, token)
                    counts[token] += 1
                else
                    counts["other"] += 1
                end
            end
            denom = max(1, nrow(raw_df))
            labels = ["none", "threads", "process", "other"]
            vals = [100.0 * counts[label] / denom for label in labels]
            p_route = Plots.bar(
                labels,
                vals;
                color=["#9aa4b2", "#3f7fb3", "#2f8f5b", "#c06c84"],
                title="Observed Outer-Route Distribution",
                xlabel="Outer Route",
                ylabel="Share of Runs [%]",
                _plot_margins(size=(2400, 760), bottom_mm=42, right_mm=42, legend=false)...
            )
            push!(panels, p_route)
        end

        if !(inner_hint_layer_df === nothing) &&
           nrow(inner_hint_layer_df) > 0 &&
           (:layer in names(inner_hint_layer_df)) &&
           (:confidence_mean in names(inner_hint_layer_df)) &&
           (:regret_mean_ns in names(inner_hint_layer_df))
            hint_df = inner_hint_layer_df[
                .!ismissing.(inner_hint_layer_df.confidence_mean) .&
                .!ismissing.(inner_hint_layer_df.regret_mean_ns), :
            ]
            if nrow(hint_df) > 0
                labels = _plot_axis_label.(String.(hint_df.layer))
                confidence = Float64.(hint_df.confidence_mean)
                regret_ms = Float64.(hint_df.regret_mean_ns) ./ 1e6
                p_hint = Plots.plot(
                    labels,
                    confidence;
                    marker=:circle,
                    color="#355070",
                    label="Mean confidence width",
                    title="Inner Adaptive Hint Confidence/Regret by Layer",
                    xlabel="Inner Layer",
                    ylabel="Confidence (lower is tighter)",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                Plots.plot!(
                    p_hint,
                    labels,
                    regret_ms;
                    marker=:diamond,
                    color="#bc4749",
                    label="Mean regret [ms]"
                )
                push!(panels, p_hint)
            end
        end

        regret_summary_df = (paper_external === nothing) ? nothing : paper_external.cross_adaptive_regret_summary_df
        if !(regret_summary_df === nothing) &&
           nrow(regret_summary_df) > 0 &&
           (:adaptive_mode in names(regret_summary_df)) &&
           (:mean_time_regret_pct in names(regret_summary_df)) &&
           (:win_rate_pct in names(regret_summary_df))
            local_df = regret_summary_df
            labels = _plot_axis_label.(String.(local_df.adaptive_mode))
            mean_regret = [v isa Missing ? NaN : Float64(v) for v in local_df.mean_time_regret_pct]
            win_rate = [v isa Missing ? NaN : Float64(v) for v in local_df.win_rate_pct]
            p_regret = Plots.plot(
                labels,
                mean_regret;
                marker=:circle,
                color="#d17a4f",
                label="Mean time regret vs best fixed [%]",
                title="Adaptive Regret vs Best Fixed (Cross-Machine)",
                xlabel="Adaptive Mode",
                ylabel="Regret / Win Rate [%]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            Plots.plot!(
                p_regret,
                labels,
                win_rate;
                marker=:utriangle,
                color="#2a9d8f",
                label="Win rate vs best fixed [%]"
            )
            push!(panels, p_regret)
        end

        if !isempty(panels)
            plt = if length(panels) == 1
                panels[1]
            else
                Plots.plot(
                    panels...;
                    layout=(length(panels), 1),
                    size=(2500, 660 * length(panels)),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
            end
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_adaptivity_behavior", spec, stamp)
        end
    end

    # 18) Paper pack: explicit layer-attribution speedup panel.
    if _paper_figure_pack_enabled() && !(paper_external === nothing)
        layer_df = paper_external.layer_attribution_speedup_df
        if !(layer_df === nothing) &&
           nrow(layer_df) > 0 &&
           (:layer_set in names(layer_df)) &&
           (:total_speedup_vs_outer_only in names(layer_df))
            speed_by_layer = Dict{String, Float64}()
            for row in eachrow(layer_df)
                layer = lowercase(strip(String(row.layer_set)))
                speed = row.total_speedup_vs_outer_only
                if speed isa Missing
                    continue
                end
                speed_f = Float64(speed)
                if !isfinite(speed_f)
                    continue
                end
                if !haskey(speed_by_layer, layer) || speed_f > speed_by_layer[layer]
                    speed_by_layer[layer] = speed_f
                end
            end
            ordered_layers = ["outer_only", "density", "thermal", "control", "multibody", "effector"]
            labels = String[]
            values = Float64[]
            for layer in ordered_layers
                if haskey(speed_by_layer, layer)
                    push!(labels, _plot_axis_label(layer))
                    push!(values, speed_by_layer[layer])
                end
            end
            if !isempty(values)
                plt = Plots.bar(
                    labels,
                    values;
                    color="#4f9d69",
                    title="Layer Attribution: Speedup vs Outer-Only Baseline",
                    xlabel="Layer Stack",
                    ylabel="Speedup vs outer-only [x]",
                    _plot_margins(size=(2300, 980), bottom_mm=60, right_mm=42, legend=false)...
                )
                Plots.hline!(plt, [1.0]; color=:black, linestyle=:dash, label=false)
                _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_layer_attribution", spec, stamp)
            end
        end
    end

    # 19) Paper pack: density-backend benchmark panel.
    if _paper_figure_pack_enabled()
        density_df = density_backend_breakdown_df === nothing ? summarize_density_backend_breakdown(raw_df) : density_backend_breakdown_df
        if nrow(density_df) > 0 &&
           (:density_backend_bucket in names(density_df)) &&
           (:total_time_mean_s in names(density_df)) &&
           (:sim_seconds_per_wall_second_mean in names(density_df))
            bucket_order = [
                "gram_point_to_point",
                "gram_surrogate",
                "gram_static_grid_or_cached_surrogate",
                "non_gram"
            ]
            labels = String[]
            mean_time = Float64[]
            throughput = Float64[]
            success_rate = Float64[]
            for bucket in bucket_order
                idx = findfirst(==(bucket), String.(density_df.density_backend_bucket))
                idx === nothing && continue
                row = density_df[idx, :]
                push!(labels, _plot_axis_label(bucket))
                push!(mean_time, row.total_time_mean_s isa Missing ? NaN : Float64(row.total_time_mean_s))
                push!(throughput, row.sim_seconds_per_wall_second_mean isa Missing ? NaN : Float64(row.sim_seconds_per_wall_second_mean))
                push!(success_rate, row.success_rate_pct isa Missing ? NaN : Float64(row.success_rate_pct))
            end
            if !isempty(labels)
                p_time = Plots.bar(
                    labels,
                    mean_time;
                    color="#7f8c8d",
                    ylabel="Mean Runtime [s]",
                    title="Density Backend Benchmark: Runtime",
                    xticks=(1:length(labels), fill("", length(labels))),
                    _plot_margins(size=(2400, 600), bottom_mm=10, right_mm=42)...
                )
                p_tp = Plots.bar(
                    labels,
                    throughput;
                    color="#2f8f5b",
                    ylabel="Sim sec / wall sec",
                    title="Density Backend Benchmark: Throughput",
                    xticks=(1:length(labels), fill("", length(labels))),
                    _plot_margins(size=(2400, 600), bottom_mm=10, right_mm=42)...
                )
                p_sr = Plots.bar(
                    labels,
                    success_rate;
                    color="#3f7fb3",
                    xlabel="Density Backend Bucket",
                    ylabel="Success Rate [%]",
                    title="Density Backend Benchmark: Solve Success",
                    _plot_margins(size=(2400, 760), bottom_mm=64, right_mm=42)...
                )
                plt = Plots.plot(
                    p_time,
                    p_tp,
                    p_sr;
                    layout=(3, 1),
                    size=(2500, 1880),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
                _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_density_backend", spec, stamp)
            end
        end
    end

    # 20) Paper pack: accuracy-parity suite.
    if _paper_figure_pack_enabled() && !(paper_external === nothing)
        panels = Any[]
        deep_df = paper_external.deep_accuracy_df
        if !(deep_df === nothing) &&
           nrow(deep_df) > 0 &&
           (:rung in names(deep_df)) &&
           (:traj_pos_rel_rms_median_pct in names(deep_df)) &&
           (:traj_vel_rel_rms_median_pct in names(deep_df))
            local_df = copy(deep_df)
            if :mode in names(local_df)
                sort!(local_df, :mode)
            end
            labels = _plot_axis_label.(String.(local_df.rung))
            pos_rms = [v isa Missing ? NaN : Float64(v) for v in local_df.traj_pos_rel_rms_median_pct]
            vel_rms = [v isa Missing ? NaN : Float64(v) for v in local_df.traj_vel_rel_rms_median_pct]
            p_traj = Plots.bar(
                labels,
                hcat(pos_rms, vel_rms);
                label=["Pos RMS rel err (median %)" "Vel RMS rel err (median %)"],
                color=["#2a9d8f" "#3f7fb3"],
                title="Accuracy Parity: Trajectory RMS Relative Error",
                xlabel="Rung",
                ylabel="Relative Error [%]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            push!(panels, p_traj)

            if (:periapsis_time_abs_err_p90_s in names(local_df)) && (:interface_time_abs_err_p90_s in names(local_df))
                peri_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.periapsis_time_abs_err_p90_s]
                interface_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.interface_time_abs_err_p90_s]
                p_event = Plots.bar(
                    labels,
                    hcat(peri_p90, interface_p90);
                    label=["Periapsis timing abs err (P90)" "Interface timing abs err (P90)"],
                    color=["#d17a4f" "#e9c46a"],
                    title="Accuracy Parity: Event-Time Error",
                    xlabel="Rung",
                    ylabel="Absolute Error [s]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_event)
            end

            if (:propellant_rel_err_p90_pct in names(local_df)) && (:control_impulse_rel_err_p90_pct in names(local_df))
                prop_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.propellant_rel_err_p90_pct]
                impulse_p90 = [v isa Missing ? NaN : Float64(v) for v in local_df.control_impulse_rel_err_p90_pct]
                p_control = Plots.bar(
                    labels,
                    hcat(prop_p90, impulse_p90);
                    label=["Propellant rel err (P90 %)" "Control impulse rel err (P90 %)"],
                    color=["#8d99ae" "#6d597a"],
                    title="Accuracy Parity: Control/Propellant",
                    xlabel="Rung",
                    ylabel="Relative Error [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_control)
            end

            if :callback_exact_match_pct in names(local_df)
                callback_match = [v isa Missing ? NaN : Float64(v) for v in local_df.callback_exact_match_pct]
                p_callback = Plots.bar(
                    labels,
                    callback_match;
                    color="#457b9d",
                    title="Accuracy Parity: Callback-State Exact Match",
                    xlabel="Rung",
                    ylabel="Exact Match [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=42, legend=false)...
                )
                push!(panels, p_callback)
            end
        end

        mc_df = paper_external.montecarlo_parity_df
        if !(mc_df === nothing) &&
           nrow(mc_df) > 0 &&
           (:mode in names(mc_df)) &&
           (:rung in names(mc_df)) &&
           (:event_time_ks_distance in names(mc_df))
            agg = combine(
                groupby(mc_df, [:mode, :rung]),
                :event_time_ks_distance => (v -> _safe_stat(v, median)) => :event_ks_median,
                :pos_ks_distance => (v -> _safe_stat(v, median)) => :pos_ks_median,
                :vel_ks_distance => (v -> _safe_stat(v, median)) => :vel_ks_median,
                :mass_ks_distance => (v -> _safe_stat(v, median)) => :mass_ks_median
            )
            sort!(agg, :mode)
            labels = _plot_axis_label.(String.(agg.rung))
            event_ks = [v isa Missing ? NaN : Float64(v) for v in agg.event_ks_median]
            pos_ks = [v isa Missing ? NaN : Float64(v) for v in agg.pos_ks_median]
            vel_ks = [v isa Missing ? NaN : Float64(v) for v in agg.vel_ks_median]
            mass_ks = [v isa Missing ? NaN : Float64(v) for v in agg.mass_ks_median]
            p_mc = Plots.bar(
                labels,
                hcat(event_ks, pos_ks, vel_ks, mass_ks);
                label=["Event-time KS" "Position KS" "Velocity KS" "Mass KS"],
                color=["#bc4749" "#2a9d8f" "#3f7fb3" "#6d597a"],
                title="Accuracy Parity: Monte Carlo Distribution KS Distance",
                xlabel="Rung",
                ylabel="KS Distance",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            push!(panels, p_mc)
        end

        if isempty(panels)
            gate_parts = DataFrame[]
            if !(split_gate_df === nothing) && nrow(split_gate_df) > 0
                local_df = copy(split_gate_df)
                local_df[!, :gate_label] = fill("split", nrow(local_df))
                push!(gate_parts, local_df)
            end
            if !(multirate_gate_df === nothing) && nrow(multirate_gate_df) > 0
                local_df = copy(multirate_gate_df)
                local_df[!, :gate_label] = fill("multirate", nrow(local_df))
                push!(gate_parts, local_df)
            end
            if !isempty(gate_parts)
                gate_df = vcat(gate_parts...; cols=:union)
                labels = _plot_axis_label.(String.(gate_df.scenario))
                pos_max = [v isa Missing ? NaN : Float64(v) for v in gate_df.pos_rel_max]
                vel_max = [v isa Missing ? NaN : Float64(v) for v in gate_df.vel_rel_max]
                p_gate = Plots.bar(
                    labels,
                    hcat(pos_max, vel_max);
                    label=["Pos rel max" "Vel rel max"],
                    color=["#2a9d8f" "#3f7fb3"],
                    title="Accuracy Gate Fallback: Trajectory Relative Error Maxima",
                    xlabel="Scenario",
                    ylabel="Relative Error",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_gate)
            end
        end

        if !isempty(panels)
            plt = if length(panels) == 1
                panels[1]
            else
                Plots.plot(
                    panels...;
                    layout=(length(panels), 1),
                    size=(2500, 660 * length(panels)),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
            end
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_accuracy_suite", spec, stamp)
        end
    end

    # 21) Paper pack: cross-machine comparison panel.
    if _paper_figure_pack_enabled() && !(paper_external === nothing)
        panels = Any[]
        speedup_summary_df = paper_external.cross_speedup_summary_df
        if !(speedup_summary_df === nothing) &&
           nrow(speedup_summary_df) > 0 &&
           (:rung in names(speedup_summary_df)) &&
           (:median_speedup_vs_r0 in names(speedup_summary_df))
            local_df = copy(speedup_summary_df)
            sort!(local_df, :median_speedup_vs_r0, rev=true)
            labels = _plot_axis_label.(String.(local_df.rung))
            med_speedup = [v isa Missing ? NaN : Float64(v) for v in local_df.median_speedup_vs_r0]
            p_speedup = Plots.bar(
                labels,
                med_speedup;
                color="#4f9d69",
                title="Cross-Machine Median Speedup vs R0",
                xlabel="Rung",
                ylabel="Median Speedup [x]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=42, legend=false)...
            )
            Plots.hline!(p_speedup, [1.0]; color=:black, linestyle=:dash, label=false)
            push!(panels, p_speedup)
        end

        regret_summary_df = paper_external.cross_adaptive_regret_summary_df
        if !(regret_summary_df === nothing) &&
           nrow(regret_summary_df) > 0 &&
           (:adaptive_mode in names(regret_summary_df)) &&
           (:mean_time_regret_pct in names(regret_summary_df))
            local_df = regret_summary_df
            labels = _plot_axis_label.(String.(local_df.adaptive_mode))
            regret_pct = [v isa Missing ? NaN : Float64(v) for v in local_df.mean_time_regret_pct]
            win_rate = (:win_rate_pct in names(local_df)) ?
                [v isa Missing ? NaN : Float64(v) for v in local_df.win_rate_pct] :
                fill(NaN, nrow(local_df))
            p_regret = Plots.plot(
                labels,
                regret_pct;
                marker=:circle,
                color="#d17a4f",
                label="Mean time regret [%]",
                title="Cross-Machine Adaptive Regret vs Best Fixed",
                xlabel="Adaptive Mode",
                ylabel="Regret / Win Rate [%]",
                _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
            )
            Plots.plot!(p_regret, labels, win_rate; marker=:diamond, color="#2a9d8f", label="Win rate [%]")
            push!(panels, p_regret)
        end

        route_mix_summary_df = paper_external.cross_route_mix_summary_df
        if !(route_mix_summary_df === nothing) &&
           nrow(route_mix_summary_df) > 0 &&
           (:none_pct_mean in names(route_mix_summary_df)) &&
           (:threads_pct_mean in names(route_mix_summary_df)) &&
           (:process_pct_mean in names(route_mix_summary_df))
            local_df = copy(route_mix_summary_df)
            if :mode in names(local_df)
                mask = [_is_adaptive_mode_token(v) for v in local_df.mode]
                if any(mask)
                    local_df = local_df[mask, :]
                end
            end
            if nrow(local_df) > 0
                labels = if :rung in names(local_df)
                    _plot_axis_label.(String.(local_df.rung))
                elseif :mode in names(local_df)
                    _plot_axis_label.(String.(local_df.mode))
                else
                    _plot_axis_label.(string.("adaptive_", collect(1:nrow(local_df))))
                end
                none_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.none_pct_mean]
                threads_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.threads_pct_mean]
                process_pct = [v isa Missing ? 0.0 : Float64(v) for v in local_df.process_pct_mean]
                p_mix = Plots.bar(
                    labels,
                    hcat(none_pct, threads_pct, process_pct);
                    label=["none" "threads" "process"],
                    bar_position=:stack,
                    color=["#9aa4b2" "#3f7fb3" "#2f8f5b"],
                    title="Cross-Machine Adaptive Route Mix",
                    xlabel="Adaptive Mode / Rung",
                    ylabel="Mean Route Share [%]",
                    _plot_margins(size=(2400, 760), bottom_mm=72, right_mm=62, legend=:outertopright)...
                )
                push!(panels, p_mix)
            end
        end

        if !isempty(panels)
            plt = if length(panels) == 1
                panels[1]
            else
                Plots.plot(
                    panels...;
                    layout=(length(panels), 1),
                    size=(2500, 660 * length(panels)),
                    left_margin=20 * Plots.mm,
                    right_margin=20 * Plots.mm
                )
            end
            _save_runtime_plot!(plot_artifacts, plt, outdir, "runtime_plot_paper_cross_machine", spec, stamp)
        end
    end

    return plot_artifacts
end

@inline function _fmt(v; digits::Int=3)
    if v isa Missing
        return "n/a"
    elseif v isa AbstractFloat
        return isfinite(v) ? string(round(v; digits=digits)) : "n/a"
    else
        return string(v)
    end
end

function _scenario_metric(summary_df::DataFrame, scenario::String, metric::Symbol)
    idx = findfirst(==(scenario), summary_df.scenario)
    idx === nothing && return nothing
    return summary_df[idx, metric]
end

function write_report(
    path::String,
    spec::ProfileSpec,
    raw_df::DataFrame,
    summary_df::DataFrame,
    orbit_summary_df::DataFrame,
    entry_duration_summary_df::Union{Nothing, DataFrame}=nothing;
    plot_paths::Vector{String}=String[],
    split_gate_df::Union{Nothing, DataFrame}=nothing,
    split_gate_csv_path::Union{Nothing, String}=nothing,
    split_gate_report_path::Union{Nothing, String}=nothing,
    multirate_gate_df::Union{Nothing, DataFrame}=nothing,
    multirate_gate_csv_path::Union{Nothing, String}=nothing,
    multirate_gate_report_path::Union{Nothing, String}=nothing,
    stage_timing_df::Union{Nothing, DataFrame}=nothing,
    inner_hint_layer_df::Union{Nothing, DataFrame}=nothing,
    inner_hint_layer_csv_path::Union{Nothing, String}=nothing,
    density_backend_breakdown_df::Union{Nothing, DataFrame}=nothing,
    density_backend_breakdown_csv_path::Union{Nothing, String}=nothing,
    entry_duration_summary_csv_path::Union{Nothing, String}=nothing
)
    generated = string(now(UTC))
    julia_ver = string(VERSION)
    nthreads = Threads.nthreads()
    solver_mode_default = _perf_default_solver_mode(spec.name)
    solver_mode_effective = _perf_solver_mode_env(spec.name)
    hw_default = _runtime_hardware_snapshot()

    function _first_nonmissing(df::DataFrame, col::Symbol, fallback)
        if !(col in names(df))
            return fallback
        end
        vals = [v for v in df[!, col] if !(v isa Missing)]
        isempty(vals) && return fallback
        return vals[1]
    end

    hardware_class = string(_first_nonmissing(raw_df, :hardware_class, hw_default.hardware_class))
    machine_label = string(_first_nonmissing(raw_df, :machine_label, hw_default.machine_label))
    host_name = string(_first_nonmissing(raw_df, :host_name, hw_default.host_name))
    cpu_name = string(_first_nonmissing(raw_df, :cpu_name, hw_default.cpu_name))
    cpu_threads = Int(_first_nonmissing(raw_df, :cpu_threads, hw_default.cpu_threads))
    julia_threads = Int(_first_nonmissing(raw_df, :julia_threads, hw_default.julia_threads))
    os_name = string(_first_nonmissing(raw_df, :os, hw_default.os))
    arch_name = string(_first_nonmissing(raw_df, :arch, hw_default.arch))

    bench_stage_s = missing
    split_stage_s = missing
    multirate_stage_s = missing
    orbit_stage_s = missing
    entry_duration_stage_s = missing
    total_stage_s = missing
    if !(stage_timing_df === nothing) && (:stage in names(stage_timing_df)) && (:elapsed_s in names(stage_timing_df))
        for row in eachrow(stage_timing_df)
            stage_name = String(row.stage)
            if stage_name == "run_benchmarks"
                bench_stage_s = row.elapsed_s
            elseif stage_name == "run_split_rollout_gate"
                split_stage_s = row.elapsed_s
            elseif stage_name == "run_multirate_rollout_gate"
                multirate_stage_s = row.elapsed_s
            elseif stage_name == "run_per_orbit"
                orbit_stage_s = row.elapsed_s
            elseif stage_name == "run_entry_duration_sweep"
                entry_duration_stage_s = row.elapsed_s
            elseif stage_name == "total"
                total_stage_s = row.elapsed_s
            end
        end
    end

    valid_rows = findall(x -> !ismissing(x), summary_df.total_time_mean_s)
    fastest = nothing
    slowest = nothing
    if !isempty(valid_rows)
        vals = summary_df.total_time_mean_s[valid_rows]
        fastest = summary_df[valid_rows[argmin(vals)], :]
        slowest = summary_df[valid_rows[argmax(vals)], :]
    end

    baseline = _scenario_metric(summary_df, PERF_BASELINE_SCENARIO, :total_time_mean_s)
    orientation = _scenario_metric(summary_df, "single_orientation_aero", :total_time_mean_s)
    multi4 = _scenario_metric(summary_df, "multi_4_gravity", :total_time_mean_s)
    multi8 = _scenario_metric(summary_df, "multi_8_gravity", :total_time_mean_s)
    multi16 = _scenario_metric(summary_df, "multi_16_gravity", :total_time_mean_s)
    multi32 = _scenario_metric(summary_df, "multi_32_gravity", :total_time_mean_s)
    multi64 = _scenario_metric(summary_df, "multi_64_gravity", :total_time_mean_s)
    harmonics20 = _scenario_metric(summary_df, "single_harmonics_l20", :total_time_mean_s)
    harmonics50 = _scenario_metric(summary_df, "single_harmonics_l50", :total_time_mean_s)

    mc_rows = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
    mc_mean = nrow(mc_rows) > 0 ? mean(mc_rows.total_time_s) : missing
    mc_p90 = nrow(mc_rows) > 0 ? quantile(mc_rows.total_time_s, 0.9) : missing
    mc_std = nrow(mc_rows) > 0 ? std(mc_rows.total_time_s; corrected=false) : missing
    fallback_rows = raw_df[(raw_df.solve_success .== true) .& (raw_df.solver_fallback_used .== true), :]
    solver_modes = _safe_unique_join(raw_df.solver_mode)
    total_success = count(identity, raw_df.solve_success)
    total_samples = nrow(raw_df)
    solve_success_rate = total_samples > 0 ? (100.0 * total_success / total_samples) : missing
    solve_failure_rate = total_samples > 0 ? (100.0 * (total_samples - total_success) / total_samples) : missing
    requested_runs = if :attempt in Symbol.(names(raw_df))
        count(v -> !(v isa Missing) && Int(v) == 1, raw_df.attempt)
    else
        total_samples
    end
    retries_total = max(total_samples - requested_runs, 0)
    retry_count_mean_all_runs = requested_runs > 0 ? (Float64(retries_total) / Float64(requested_runs)) : missing
    total_attempt_time_s = begin
        vals = [Float64(v) for v in raw_df.total_time_s if !(v isa Missing)]
        isempty(vals) ? missing : sum(vals)
    end
    penalized_expected_wall_time_all_runs = (total_attempt_time_s isa Missing || requested_runs <= 0) ? missing :
        (Float64(total_attempt_time_s) / Float64(requested_runs))
    fallback_count_mean_all_attempts_global = if :solver_fallback_count in Symbol.(names(raw_df))
        vals = Float64[(v isa Missing) ? 0.0 : Float64(v) for v in raw_df.solver_fallback_count]
        isempty(vals) ? missing : mean(vals)
    else
        missing
    end
    robustness_table = DataFrame()
    if (:penalized_expected_wall_time_s in Symbol.(names(summary_df))) &&
       (:failure_rate in Symbol.(names(summary_df))) &&
       (:retry_count_mean in Symbol.(names(summary_df)))
        robustness_table = summary_df[.!ismissing.(summary_df.penalized_expected_wall_time_s), :]
        if nrow(robustness_table) > 0
            sort!(robustness_table, :penalized_expected_wall_time_s)
        end
    end
    failed_groups = summary_df[summary_df.samples_failed .> 0, :]

    spice_rows = summary_df[.!ismissing.(summary_df.spice_calls_total_mean), :]
    spice_peak = nothing
    if nrow(spice_rows) > 0
        spice_rows_sorted = copy(spice_rows)
        sort!(spice_rows_sorted, :spice_calls_total_mean, rev=true)
        spice_peak = spice_rows_sorted[1, :]
    end

    ci_rows = summary_df[
        .!ismissing.(summary_df.total_time_ci95_low_s) .&
        .!ismissing.(summary_df.total_time_ci95_high_s) .&
        .!ismissing.(summary_df.total_time_mean_s), :
    ]
    ci_rows_sorted = DataFrame()
    if nrow(ci_rows) > 0
        ci_rows_sorted = copy(ci_rows)
        ci_rows_sorted[!, :_total_sort] = [Float64(v) for v in ci_rows_sorted.total_time_mean_s]
        sort!(ci_rows_sorted, :_total_sort, rev=true)
        select!(ci_rows_sorted, Not(:_total_sort))
    end

    split_pass_count = 0
    split_total = 0
    split_any_fail = false
    if !(split_gate_df === nothing) && (:pass_all in names(split_gate_df))
        split_total = nrow(split_gate_df)
        split_pass_count = count(Bool.(split_gate_df.pass_all))
        split_any_fail = split_pass_count < split_total
    end

    multirate_pass_count = 0
    multirate_total = 0
    multirate_any_fail = false
    if !(multirate_gate_df === nothing) && (:pass_all in names(multirate_gate_df))
        multirate_total = nrow(multirate_gate_df)
        multirate_pass_count = count(Bool.(multirate_gate_df.pass_all))
        multirate_any_fail = multirate_pass_count < multirate_total
    end

    open(path, "w") do io
        println(io, "# SpaceAGORA Computational Time Analysis (`$(spec.name)` profile)")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$julia_ver`")
        println(io, "- Threads: `$nthreads`")
        println(io, "- Machine label: `$(machine_label)`")
        println(io, "- Hardware class: `$(hardware_class)`")
        println(io, "- Hostname: `$(host_name)`")
        println(io, "- CPU: `$(cpu_name)` (`$(cpu_threads)` system threads)")
        println(io, "- Julia threads in process: `$(julia_threads)`")
        println(io, "- OS/Arch: `$(os_name)` / `$(arch_name)`")
        println(io, "- Repeats per deterministic scenario: `$(spec.repeats)`")
        println(io, "- Warmup runs per scenario: `$(spec.warmup)`")
        println(io, "- Monte Carlo seeds: `$(spec.montecarlo_samples)`")
        println(io, "- Solver mode default for profile: `$(solver_mode_default)`")
        println(io, "- Solver mode effective (after env overrides): `$(solver_mode_effective)`")
        println(io, "- Solver mode(s) observed: `$(_fmt(solver_modes))`")
        if !(total_stage_s isa Missing)
            println(
                io,
                "- Stage elapsed [s]: run_benchmarks=`$(_fmt(bench_stage_s))`, " *
                "split_gate=`$(_fmt(split_stage_s))`, multirate_gate=`$(_fmt(multirate_stage_s))`, " *
                "mission_time_sweep=`$(_fmt(orbit_stage_s))`, entry_duration_sweep=`$(_fmt(entry_duration_stage_s))`, " *
                "total=`$(_fmt(total_stage_s))`"
            )
        end
        println(io)
        println(io, "## Claim Scope")
        println(io)
        if spec.name == "full"
            println(io, "- This `full` profile report is intended for paper-grade scalability and adaptivity claims.")
        else
            println(io, "- This `quick` profile report is for development/CI/regression and should not be used for paper-grade scalability claims.")
        end
        println(io, "- Mission-time sweep excludes `entry` cases by design because multipliers are based on baseline orbital periods.")
        println(io, "- Entry behavior is covered separately by the `Entry-Duration Sweep` section using atmospheric-interface counts.")
        println(io)
        println(io, "## Key Findings")
        println(io)
        if fastest === nothing || slowest === nothing
            println(io, "- No successful runs were recorded.")
        else
            println(io, "- Slowest successful scenario: `$(slowest.scenario)` with mean total time `$(round(slowest.total_time_mean_s; digits=3)) s`.")
            println(io, "- Fastest successful scenario: `$(fastest.scenario)` with mean total time `$(round(fastest.total_time_mean_s; digits=3)) s`.")
        end
        if baseline !== nothing && orientation !== nothing && !ismissing(baseline) && !ismissing(orientation) && baseline > 0.0
            println(io, "- Orientation + aerodynamic run vs `$(PERF_BASELINE_SCENARIO)`: `$(round(orientation / baseline; digits=2))x` runtime.")
        end
        if multi4 !== nothing && multi8 !== nothing && multi16 !== nothing && multi32 !== nothing && multi64 !== nothing &&
           !ismissing(multi4) && !ismissing(multi8) && !ismissing(multi16) && !ismissing(multi32) && !ismissing(multi64) && multi4 > 0.0
            println(
                io,
                "- Multi-satellite scaling (runtime): " *
                "`8/4=$(round(multi8 / multi4; digits=2))x`, " *
                "`16/4=$(round(multi16 / multi4; digits=2))x`, " *
                "`32/4=$(round(multi32 / multi4; digits=2))x`, " *
                "`64/4=$(round(multi64 / multi4; digits=2))x`."
            )
        end
        if harmonics20 !== nothing && harmonics50 !== nothing && !ismissing(harmonics20) && !ismissing(harmonics50) && harmonics20 > 0.0
            println(io, "- Harmonics scaling: `L=50` is `$(round(harmonics50 / harmonics20; digits=2))x` relative to `L=20`.")
        end
        if !(mc_mean isa Missing)
            println(io, "- Monte Carlo runtime spread: mean `$(round(mc_mean; digits=3)) s`, p90 `$(round(mc_p90; digits=3)) s`, std `$(round(mc_std; digits=3)) s`.")
        end
        println(io, "- Auto-stiff fallback activations (successful runs): `$(nrow(fallback_rows))`.")
        println(io, "- Successful samples: `$(total_success)/$(total_samples)` (`$(_fmt(solve_success_rate))%`).")
        println(io, "- Failed attempts: `$(total_samples - total_success)/$(total_samples)` (`$(_fmt(solve_failure_rate))%`).")
        println(io, "- Retry overhead: `$(retries_total)` retries across `$(requested_runs)` requested runs (`$(_fmt(retry_count_mean_all_runs))` retries/requested run).")
        println(io, "- Robustness-adjusted expected wall time (all attempts): `$(_fmt(penalized_expected_wall_time_all_runs)) s/requested run`.")
        println(io, "- Mean fallback count across all attempts: `$(_fmt(fallback_count_mean_all_attempts_global))`.")
        if !(spice_peak === nothing)
            println(
                io,
                "- Highest mean SPICE call budget: `$(spice_peak.scenario)` with `$(_fmt(spice_peak.spice_calls_total_mean))` calls/run " *
                "(`$(_fmt(spice_peak.spice_calls_per_wall_second_mean))` calls/wall-sec)."
            )
        end
        if split_total > 0
            println(io, "- Split rollout guardrail: `$(split_pass_count)/$(split_total)` pass.")
        else
            println(io, "- Split rollout guardrail: disabled or no gate rows.")
        end
        if multirate_total > 0
            println(io, "- Multirate rollout guardrail: `$(multirate_pass_count)/$(multirate_total)` pass.")
        else
            println(io, "- Multirate rollout guardrail: disabled or no gate rows.")
        end
        println(io, "- Plot artifacts generated: `$(length(plot_paths))`.")
        if nrow(failed_groups) > 0
            println(io, "- Solver failures detected in `$(nrow(failed_groups))` scenario groups; timings only use successful runs.")
        end
        if !isempty(plot_paths)
            println(io)
            println(io, "## Plot Artifacts")
            println(io)
            for plot_path in plot_paths
                println(io, "- `$(plot_path)`")
            end
            println(io)
        end

        println(io, "## Statistical Confidence")
        println(io)
        if nrow(ci_rows_sorted) == 0
            println(io, "- No confidence interval rows are available (insufficient successful samples).")
        else
            println(io, "| Scenario | Samples | Mean Total (s) | 95% CI Low (s) | 95% CI High (s) | SEM (s) | CV (%) |")
            println(io, "|---|---:|---:|---:|---:|---:|---:|")
            for row in eachrow(ci_rows_sorted)
                println(
                    io,
                    "| $(row.scenario) | $(row.samples_success) | $(_fmt(row.total_time_mean_s)) | $(_fmt(row.total_time_ci95_low_s)) | " *
                    "$(_fmt(row.total_time_ci95_high_s)) | $(_fmt(row.total_time_sem_s)) | $(_fmt(row.total_time_cv_pct)) |"
                )
            end
        end
        println(io)

        println(io, "## Subsystem Evidence (SPICE)")
        println(io)
        if nrow(spice_rows) == 0
            println(io, "- No SPICE counters were captured for successful rows.")
        else
            println(io, "| Scenario | N-body total calls | SRP total calls | Planet pxform calls | SPICE total calls | Calls/wall-sec |")
            println(io, "|---|---:|---:|---:|---:|---:|")
            spice_table = copy(spice_rows)
            sort!(spice_table, :spice_calls_total_mean, rev=true)
            for row in eachrow(spice_table)
                println(
                    io,
                    "| $(row.scenario) | $(_fmt(row.nbody_spkpos_total_calls_mean)) | $(_fmt(row.srp_spkpos_total_calls_mean)) | " *
                    "$(_fmt(row.planet_pxform_total_calls_mean)) | $(_fmt(row.spice_calls_total_mean)) | $(_fmt(row.spice_calls_per_wall_second_mean)) |"
                )
            end
        end
        println(io)

        println(io, "## Inner Hint Layer Stats")
        println(io)
        hint_rows = inner_hint_layer_df === nothing ? 0 : nrow(inner_hint_layer_df)
        if hint_rows == 0
            println(io, "- No persistent inner-hint rows matched the active parallel profile/machine filter.")
        else
            println(io, "| Layer | Signatures | Choices | Samples | Mean elapsed (ns) | Mean confidence | Mean regret (ns) |")
            println(io, "|---|---:|---:|---:|---:|---:|---:|")
            hint_table = copy(inner_hint_layer_df)
            sort!(hint_table, :regret_mean_ns, rev=true)
            for row in eachrow(hint_table)
                println(
                    io,
                    "| $(row.layer) | $(row.signature_count) | $(row.choice_count) | $(row.samples_total) | " *
                    "$(_fmt(row.elapsed_mean_ns)) | $(_fmt(row.confidence_mean)) | $(_fmt(row.regret_mean_ns)) |"
                )
            end
        end
        if !(inner_hint_layer_csv_path === nothing)
            println(io, "- Inner hint CSV: `$(inner_hint_layer_csv_path)`")
        end
        println(io)

        println(io, "## Density Backend Parallelism Scope")
        println(io)
        println(io, "- Callback-level parallelism is available for density callbacks.")
        println(io, "- True density-backend scalability depends on the active density model path.")
        println(io, "- GRAM point-to-point is lock-sensitive and generally prefers outer process isolation for heavy workloads.")
        println(io, "- For throughput-heavy campaigns, prefer surrogate/static-grid/cached-surrogate paths (batched/vectorized surrogates where available).")
        density_rows = density_backend_breakdown_df === nothing ? 0 : nrow(density_backend_breakdown_df)
        if density_rows == 0
            println(io, "- No density-backend breakdown rows were produced.")
        else
            println(io, "| Density Backend Bucket | Covered | Density Families | Outer Routes | Success/Total | Success Rate (%) | Mean Total (s) | P90 Total (s) | Mean Sim sec / wall sec | Recommended Route |")
            println(io, "|---|---|---|---|---:|---:|---:|---:|---:|---|")
            for row in eachrow(density_backend_breakdown_df)
                println(
                    io,
                    "| $(row.density_backend_bucket) | $(row.covered) | $(_fmt(row.density_families)) | $(_fmt(row.outer_routes)) | " *
                    "$(row.samples_success)/$(row.samples_total) | $(_fmt(row.success_rate_pct)) | " *
                    "$(_fmt(row.total_time_mean_s)) | $(_fmt(row.total_time_p90_s)) | " *
                    "$(_fmt(row.sim_seconds_per_wall_second_mean)) | $(row.recommended_route) |"
                )
            end
        end
        if !(density_backend_breakdown_csv_path === nothing)
            println(io, "- Density backend breakdown CSV: `$(density_backend_breakdown_csv_path)`")
        end
        println(io)

        println(io, "## Fidelity Guardrails")
        println(io)
        println(io, "- Split rollout gate enabled: `$(_split_rollout_enabled())`")
        println(io, "- Split rollout enforce mode: `$(_split_rollout_enforce())`")
        println(io, "- Split rollout cases: `$(join(_split_rollout_case_names(), ", "))`")
        println(io, "- Split rollout solvers: `$(join(_split_rollout_solver_variants(), ", "))`")
        if split_total > 0
            println(io, "- Gate pass count: `$(split_pass_count)/$(split_total)`")
            println(io, "- Any gate failure: `$(split_any_fail)`")
            if !(split_gate_csv_path === nothing)
                println(io, "- Gate CSV: `$(split_gate_csv_path)`")
            end
            if !(split_gate_report_path === nothing)
                println(io, "- Gate report: `$(split_gate_report_path)`")
            end
        else
            println(io, "- Gate rows: none (gate disabled or no matching cases).")
        end
        println(io)
        println(io, "- Multirate rollout gate enabled: `$(_multirate_rollout_enabled())`")
        println(io, "- Multirate rollout enforce mode: `$(_multirate_rollout_enforce())`")
        println(io, "- Multirate rollout cases: `$(join(_multirate_rollout_case_names(), ", "))`")
        println(io, "- Multirate rollout max slowdown ratio: `$(_multirate_rollout_max_slowdown_ratio())`")
        if multirate_total > 0
            println(io, "- Multirate gate pass count: `$(multirate_pass_count)/$(multirate_total)`")
            println(io, "- Multirate any gate failure: `$(multirate_any_fail)`")
            if !(multirate_gate_csv_path === nothing)
                println(io, "- Multirate gate CSV: `$(multirate_gate_csv_path)`")
            end
            if !(multirate_gate_report_path === nothing)
                println(io, "- Multirate gate report: `$(multirate_gate_report_path)`")
            end
        else
            println(io, "- Multirate gate rows: none (gate disabled or no matching cases).")
        end
        println(io)

        println(io, "## Verification Snapshot")
        println(io)
        println(io, "- Solver-success samples: `$(total_success)/$(total_samples)` (`$(_fmt(solve_success_rate))%`).")
        println(io, "- Scenario groups with failures: `$(nrow(failed_groups))`.")
        println(io, "- Auto-stiff fallback rows: `$(nrow(fallback_rows))`.")
        if split_total > 0
            println(io, "- Split rollout verification rows: `$(split_total)` (`$(split_pass_count)` pass).")
        end
        if multirate_total > 0
            println(io, "- Multirate rollout verification rows: `$(multirate_total)` (`$(multirate_pass_count)` pass).")
        end
        println(io)

        println(io)
        println(io, "## Robustness-Adjusted Scenario Summary")
        println(io)
        if nrow(robustness_table) == 0
            println(io, "- No robustness-adjusted rows are available.")
        else
            println(io, "| Scenario | Success/Total | Requested Runs | Retries Total | Retry Mean | Failure Rate (%) | Fallback Mean (All Attempts) | Penalized Expected Wall Time (s/requested run) |")
            println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
            for row in eachrow(robustness_table)
                failure_pct = row.failure_rate isa Missing ? missing : (100.0 * Float64(row.failure_rate))
                println(
                    io,
                    "| $(row.scenario) | $(row.samples_success)/$(row.samples_total) | $(row.requested_runs) | " *
                    "$(row.retries_total) | $(_fmt(row.retry_count_mean)) | $(_fmt(failure_pct)) | " *
                    "$(_fmt(row.fallback_count_mean_all_attempts)) | $(_fmt(row.penalized_expected_wall_time_s)) |"
                )
            end
        end
        println(io)

        println(io)
        println(io, "## Scenario Summary")
        println(io)
        println(io, "| Scenario | Category | Success/Total | Solver(s) | Fallback Any | Mean Total (s) | 95% CI [low, high] (s) | CV (%) | P90 (s) | Mean Solve (s) | Mean Copy (s) | SPICE Calls/run | Sim sec / wall sec | Rel. Baseline |")
        println(io, "|---|---|---:|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in eachrow(summary_df)
            ci_band = "[$(_fmt(row.total_time_ci95_low_s)), $(_fmt(row.total_time_ci95_high_s))]"
            println(
                io,
                "| $(row.scenario) | $(row.category) | $(row.samples_success)/$(row.samples_total) | $(_fmt(row.solver_sequences)) | $(_fmt(row.solver_fallback_any)) | $(_fmt(row.total_time_mean_s)) | $(ci_band) | $(_fmt(row.total_time_cv_pct)) | $(_fmt(row.total_time_p90_s)) | $(_fmt(row.solve_time_mean_s)) | $(_fmt(row.copy_time_mean_s)) | $(_fmt(row.spice_calls_total_mean)) | $(_fmt(row.sim_seconds_per_wall_second_mean)) | $(_fmt(row.relative_to_baseline)) |"
            )
        end
        println(io)
        println(io, "## Mission-Time Sweep Results (All Scenarios)")
        println(io)
        println(io, "_Entry scenarios are intentionally excluded from this sweep because multipliers are baseline-orbit based._")
        println(io)
        println(io, "| Scenario | Category | Mission-Time Multiplier (x baseline period) | Success/Total | Mission Time (s) | Mean Total (s) | 95% CI [low, high] (s) | P90 (s) | Time / Baseline-Period Unit (s) | Baseline-Period Units / Wall-sec |")
        println(io, "|---|---|---:|---:|---:|---:|---|---:|---:|---:|")
        for row in eachrow(orbit_summary_df)
            ci_band = "[$(_fmt(row.total_time_ci95_low_s)), $(_fmt(row.total_time_ci95_high_s))]"
            multiplier = hasproperty(row, :mission_time_multiplier) ? row.mission_time_multiplier : row.orbit_count
            time_per_unit = hasproperty(row, :time_per_baseline_period_mean_s) ? row.time_per_baseline_period_mean_s : row.time_per_orbit_mean_s
            units_per_wall = hasproperty(row, :baseline_periods_per_wall_second_mean) ? row.baseline_periods_per_wall_second_mean : row.orbits_per_wall_second_mean
            println(
                io,
                "| $(row.scenario) | $(row.category) | $(multiplier) | $(row.samples_success)/$(row.samples_total) | $(_fmt(row.mission_time_mean_s)) | $(_fmt(row.total_time_mean_s)) | $(ci_band) | $(_fmt(row.total_time_p90_s)) | $(_fmt(time_per_unit)) | $(_fmt(units_per_wall)) |"
            )
        end
        println(io)
        println(io, "## Entry-Duration Sweep Results")
        println(io)
        if entry_duration_summary_df === nothing || nrow(entry_duration_summary_df) == 0
            println(io, "- No entry-duration sweep rows were produced.")
        else
            measured_df = if :entry_run_role in names(entry_duration_summary_df)
                entry_duration_summary_df[entry_duration_summary_df.entry_run_role .== "measured", :]
            else
                entry_duration_summary_df
            end
            if nrow(measured_df) == 0
                println(io, "- Entry-duration sweep executed, but no measured rows were available.")
            else
                println(io, "| Scenario | Atmos Interfaces | Success/Total | Passage Duration Mean (s) | Event-Time Abs Error Mean (s) | Wall-Time / Passage Mean (s) | Mean Total (s) |")
                println(io, "|---|---:|---:|---:|---:|---:|---:|")
                for row in eachrow(measured_df)
                    println(
                        io,
                        "| $(row.scenario) | $(row.entry_atmospheric_interface_count) | $(row.samples_success)/$(row.samples_total) | " *
                        "$(_fmt(row.passage_duration_mean_s)) | $(_fmt(row.event_time_abs_error_mean_s)) | " *
                        "$(_fmt(row.wall_time_per_passage_mean_s)) | $(_fmt(row.total_time_mean_s)) |"
                    )
                end
            end
            if !(entry_duration_summary_csv_path === nothing)
                println(io)
                println(io, "- Entry-duration summary CSV: `$(entry_duration_summary_csv_path)`")
            end
        end
    end
end

