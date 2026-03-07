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

