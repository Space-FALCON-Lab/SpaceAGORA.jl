using Dates
using CSV
using DataFrames
using Statistics

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: REPO_ROOT, env_path, resolve_output_path, stored_output_path

const SIM_FILE = joinpath(@__DIR__, "..", "..", "navigation.jl")
const MC_OUTPUT_ROOT = env_path(
    "SPACEAGORA_MC_OUTPUT",
    "nominal"
)
const MC_RUNS = parse(Int, get(ENV, "SPACEAGORA_MC_RUNS", "100"))
const MC_FIRST_RUN = parse(Int, get(ENV, "SPACEAGORA_MC_FIRST_RUN", "1"))
const MC_BASE_SEED = parse(Int, get(ENV, "SPACEAGORA_MC_BASE_SEED", "202607180"))
const MC_SEED_STRIDE = parse(Int, get(ENV, "SPACEAGORA_MC_SEED_STRIDE", "10"))
const MC_RESUME = lowercase(strip(get(ENV, "SPACEAGORA_MC_RESUME", "true"))) in
    ("1", "true", "yes", "y", "on")
const MC_AGGREGATE_EVERY = parse(Int, get(ENV, "SPACEAGORA_MC_AGGREGATE_EVERY", "10"))
const MC_METRICS_SCHEMA_VERSION = 5
const DEFAULT_MISSION_TIME_SEC = "10000.0"
const DEFAULT_TARGET_COUNT = "300"
const DEFAULT_HEAP_SIZE_HINT = "4G"

@inline _truthy(value::AbstractString) = lowercase(strip(value)) in ("1", "true", "yes", "y", "on")

# When enabled, new runs log exact truth geometry at every IOD event and the
# campaign creates <output>/iod_geometry_analysis/*.csv. Legacy campaigns are
# post-processed by a truth-only replay with the full simulation dynamics.
const MC_IOD_GEOMETRY_ANALYSIS = _truthy(get(
    ENV,
    "SPACEAGORA_MC_IOD_GEOMETRY_ANALYSIS",
    "false"
))
const MC_SAVE_IOD_GEOMETRY = _truthy(get(
    ENV,
    "SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY",
    string(MC_IOD_GEOMETRY_ANALYSIS)
))
const MC_SAVE_IOD_PAIRWISE = _truthy(get(
    ENV,
    "SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS",
    "false"
))
const MC_IOD_GEOMETRY_REPLAY = _truthy(get(
    ENV,
    "SPACEAGORA_MC_IOD_GEOMETRY_REPLAY",
    "true"
))
const MC_IOD_GEOMETRY_MAX_RUNS = parse(
    Int,
    get(ENV, "SPACEAGORA_MC_IOD_GEOMETRY_MAX_RUNS", "0")
)
MC_IOD_GEOMETRY_ANALYSIS && include(joinpath(
    @__DIR__,
    "..",
    "analysis",
    "iod_geometry.jl"
))

function _selected_cases()::Vector{String}
    default_cases = join((
        "proposed",
        "centralized_oracle",
        "independent_local_da",
        "distributed_oracle_da"
    ), ',')
    raw = strip(get(ENV, "SPACEAGORA_MC_CASES", default_cases))
    cases = [strip(value) for value in split(raw, ',') if !isempty(strip(value))]
    isempty(cases) && error("SPACEAGORA_MC_CASES must contain at least one case")
    valid = Set((
        "proposed", "no_da", "centralized_oracle", "independent_local_da",
        "distributed_oracle_da", "baseline_da",
        "baseline_da_theta_0p01", "baseline_da_theta_0p05"
    ))
    invalid = [case for case in cases if !(case in valid)]
    isempty(invalid) || error("Invalid SPACEAGORA_MC_CASES entries: $(invalid)")
    return unique(cases)
end

@inline function _nav_case(case::String)::String
    startswith(case, "baseline_da_theta_") && return "baseline_da"
    return case
end

@inline function _baseline_theta(case::String)::Union{Nothing, Float64}
    case == "baseline_da_theta_0p01" && return 0.01
    case == "baseline_da_theta_0p05" && return 0.05
    return nothing
end

@inline function _run_seeds(realization::Int)
    scenario = MC_BASE_SEED + MC_SEED_STRIDE * (realization - 1)
    return (scenario=scenario, sensor=scenario + 1, od=scenario + 2, bias=scenario + 3)
end

@inline _run_dir(realization::Int) = joinpath(MC_OUTPUT_ROOT, "run_$(lpad(realization, 4, '0'))")
@inline _case_dir(realization::Int, case::String) = joinpath(_run_dir(realization), case)
@inline _metrics_path(realization::Int, case::String) = joinpath(_case_dir(realization, case), "association_quality_table.csv")

function _completed_with_seed(realization::Int, case::String, scenario_seed::Int)::Bool
    path = _metrics_path(realization, case)
    isfile(path) || return false
    try
        table = CSV.read(path, DataFrame)
        rows = table[(table.section .== "run_config") .& (table.metric .== "scenario_seed"), :]
        schema_rows = table[
            (table.section .== "run_config") .& (table.metric .== "mc_metrics_schema_version"),
            :
        ]
        complete = nrow(rows) == 1 &&
            round(Int, rows.value[1]) == scenario_seed &&
            nrow(schema_rows) == 1 &&
            round(Int, schema_rows.value[1]) == MC_METRICS_SCHEMA_VERSION
        complete || return false
        if MC_SAVE_IOD_PAIRWISE
            isfile(joinpath(
                _case_dir(realization, case),
                "iod_pairwise_consistency_table.csv"
            )) || return false
        end
        return true
    catch
        return false
    end
end

function _child_environment(realization::Int, case::String, seeds)
    child_env = copy(ENV)
    child_env["SPACEAGORA_COMPARISON_OUTPUT"] = _run_dir(realization)
    child_env["SPACEAGORA_NAV_CASE"] = _nav_case(case)
    child_env["SPACEAGORA_NAV_OUTPUT_LABEL"] = case
    baseline_theta = _baseline_theta(case)
    if baseline_theta !== nothing
        child_env["SPACEAGORA_BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD"] = string(baseline_theta)
    end
    child_env["SPACEAGORA_SCENARIO_SEED"] = string(seeds.scenario)
    child_env["SPACEAGORA_SENSOR_SEED"] = string(seeds.sensor)
    child_env["SPACEAGORA_OBSERVER_OD_SEED"] = string(seeds.od)
    child_env["SPACEAGORA_MEASUREMENT_BIAS_SEED"] = string(seeds.bias)
    child_env["SPACEAGORA_MISSION_TIME_SEC"] = get(ENV, "SPACEAGORA_MISSION_TIME_SEC", DEFAULT_MISSION_TIME_SEC)
    child_env["SPACEAGORA_N_CLUSTER_TARGETS"] = get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", DEFAULT_TARGET_COUNT)
    child_env["SPACEAGORA_SAVE_SIMULATION_RESULTS"] = "false"
    child_env["SPACEAGORA_SAVE_TARGET_ESTIMATE_FIELDS"] = "false"
    child_env["SPACEAGORA_SAVE_BUNDLE"] = "0"
    child_env["SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES"] = get(
        ENV,
        "SPACEAGORA_MC_SAVE_DETAILED_TABLES",
        "false"
    )
    child_env["SPACEAGORA_ENABLE_NAV_TIMING"] = get(ENV, "SPACEAGORA_MC_ENABLE_NAV_TIMING", "false")
    child_env["SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY"] =
        string(MC_SAVE_IOD_GEOMETRY)
    child_env["SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS"] =
        string(MC_SAVE_IOD_PAIRWISE)
    return child_env
end

function _failure_excerpt(log_path::String)::String
    isfile(log_path) || return "run log was not created"
    lines = readlines(log_path)
    isempty(lines) && return "run log is empty"
    error_index = findlast(line ->
        occursin("ERROR:", line) || occursin("LoadError", line) ||
        occursin("ArgumentError", line) || occursin("UndefVarError", line),
        lines
    )
    if error_index !== nothing
        return strip(join(lines[error_index:min(error_index + 2, length(lines))], " | "))
    end
    first_nonempty = findfirst(line -> !isempty(strip(line)), Iterators.reverse(lines))
    first_nonempty === nothing && return "run log contains only blank lines"
    tail = [strip(line) for line in lines[max(1, length(lines) - 4):end] if !isempty(strip(line))]
    return join(tail, " | ")
end

function _run_one(realization::Int, case::String)
    seeds = _run_seeds(realization)
    case_dir = _case_dir(realization, case)
    mkpath(case_dir)
    log_path = joinpath(case_dir, "run.log")

    if MC_RESUME && _completed_with_seed(realization, case, seeds.scenario)
        return (
            realization=realization, case=case, status="resumed", exit_code=0,
            elapsed_sec=0.0, scenario_seed=seeds.scenario, sensor_seed=seeds.sensor,
            observer_od_seed=seeds.od, measurement_bias_seed=seeds.bias,
            results_dir=stored_output_path(case_dir, MC_OUTPUT_ROOT),
            log_path=stored_output_path(log_path, MC_OUTPUT_ROOT)
        )
    end

    start_time = time()
    child_env = _child_environment(realization, case, seeds)
    julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
    heap_hint = strip(get(ENV, "SPACEAGORA_JULIA_HEAP_SIZE_HINT", DEFAULT_HEAP_SIZE_HINT))
    cmd = isempty(heap_hint) ?
        setenv(`$julia_exe --project=$REPO_ROOT $SIM_FILE`, child_env) :
        setenv(`$julia_exe --heap-size-hint=$heap_hint --project=$REPO_ROOT $SIM_FILE`, child_env)

    println("[$(realization)/$(MC_FIRST_RUN + MC_RUNS - 1)] $(case): seeds $(seeds)")
    flush(stdout)
    exit_code = 1
    open(log_path, "w") do log_io
        process = run(pipeline(ignorestatus(cmd); stdout=log_io, stderr=log_io); wait=false)
        last_heartbeat = time()
        while process_running(process)
            sleep(1.0)
            if time() - last_heartbeat >= 30.0
                elapsed = time() - start_time
                println(
                    "  realization=$(realization), case=$(case): still running " *
                    "($(round(elapsed / 60; digits=1)) min), log=$(log_path)"
                )
                flush(stdout)
                last_heartbeat = time()
            end
        end
        wait(process)
        exit_code = try Int(process.exitcode) catch; success(process) ? 0 : 1 end
    end
    ok = exit_code == 0 && isfile(_metrics_path(realization, case))
    if ok
        println("  completed in $(round(time() - start_time; digits=1)) s")
    else
        println("  FAILED (exit_code=$(exit_code)): $(_failure_excerpt(log_path))")
    end
    flush(stdout)
    return (
        realization=realization, case=case, status=ok ? "ok" : "failed",
        exit_code=exit_code, elapsed_sec=time() - start_time,
        scenario_seed=seeds.scenario, sensor_seed=seeds.sensor,
        observer_od_seed=seeds.od, measurement_bias_seed=seeds.bias,
        results_dir=stored_output_path(case_dir, MC_OUTPUT_ROOT),
        log_path=stored_output_path(log_path, MC_OUTPUT_ROOT)
    )
end

function _upsert_status!(status::DataFrame, row)::DataFrame
    ncol(status) == 0 && return DataFrame([row])
    matching_rows = findall(
        (status.realization .== row.realization) .&
        (status.case .== row.case)
    )
    if isempty(matching_rows)
        push!(status, row)
    else
        status[first(matching_rows), :] = row
        length(matching_rows) > 1 && deleteat!(status, matching_rows[2:end])
    end
    return status
end

function _read_all_metrics(status::DataFrame)::DataFrame
    tables = DataFrame[]
    for row in eachrow(status)
        row.status in ("ok", "resumed") || continue
        results_dir = resolve_output_path(String(row.results_dir), MC_OUTPUT_ROOT)
        path = joinpath(results_dir, "association_quality_table.csv")
        isfile(path) || continue
        table = CSV.read(path, DataFrame)
        insertcols!(table, 1,
            :realization => fill(Int(row.realization), nrow(table)),
            :case => fill(String(row.case), nrow(table)),
            :scenario_seed => fill(Int(row.scenario_seed), nrow(table)),
            :sensor_seed => fill(Int(row.sensor_seed), nrow(table)),
            :observer_od_seed => fill(Int(row.observer_od_seed), nrow(table))
        )
        push!(tables, table)
    end
    return isempty(tables) ? DataFrame() : vcat(tables...; cols=:union)
end

function _combine_case_table(status::DataFrame, filename::String; one_case_per_realization::Bool=false)::DataFrame
    tables = DataFrame[]
    seen_realizations = Set{Int}()
    for row in eachrow(status)
        row.status in ("ok", "resumed") || continue
        realization = Int(row.realization)
        one_case_per_realization && realization in seen_realizations && continue
        results_dir = resolve_output_path(String(row.results_dir), MC_OUTPUT_ROOT)
        path = joinpath(results_dir, filename)
        isfile(path) || continue
        table = CSV.read(path, DataFrame)
        insertcols!(table, 1,
            :realization => fill(realization, nrow(table)),
            :case => fill(String(row.case), nrow(table)),
            :scenario_seed => fill(Int(row.scenario_seed), nrow(table))
        )
        push!(tables, table)
        one_case_per_realization && push!(seen_realizations, realization)
    end
    return isempty(tables) ? DataFrame() : vcat(tables...; cols=:union)
end

function _metric_summary(long::DataFrame)::DataFrame
    rows = NamedTuple[]
    nrow(long) == 0 && return DataFrame(rows)
    for group in groupby(long, [:case, :section, :metric])
        values = [Float64(value) for value in group.value if isfinite(Float64(value))]
        isempty(values) && continue
        push!(rows, (
            case=String(group.case[1]), section=String(group.section[1]), metric=String(group.metric[1]),
            count=length(values), mean=mean(values), std=length(values) > 1 ? std(values) : NaN,
            minimum=minimum(values), p05=quantile(values, 0.05), median=median(values),
            p95=quantile(values, 0.95), maximum=maximum(values)
        ))
    end
    return DataFrame(rows)
end

function _problem_characterization(wide::DataFrame)::DataFrame
    nrow(wide) == 0 && return DataFrame()
    wanted = [
        "association_health.identity_anomaly_total",
        "meas_assoc.commit_accuracy_pct", "meas_assoc.recall_pct",
        "track_assoc.tt_accuracy_pct_known_only", "track_assoc.tt_recall_pct",
        "cross_m2m.iod_groups_mixed_target", "track_assoc.consensus_group_mixed_target",
        "track_lifecycle.id_switch_total", "track_lifecycle.fragment_excess_total",
        "track_lifecycle.initialization_success_pct",
        "track_lifecycle.never_correctly_initialized_unique_targets",
        "tracking.success_rate_possible_pct", "tracking.converged_mean_error_m",
        "track_lifecycle.initialization_position_error_mean_m",
        "geometry.max_simultaneously_visible_targets", "geometry.minimum_angular_separation_deg"
    ]
    columns = [:realization, :case, :scenario_seed, :sensor_seed, :observer_od_seed]
    append!(columns, [Symbol(name) for name in wanted if Symbol(name) in propertynames(wide)])
    result = select(wide, columns)
    anomaly = Symbol("association_health.identity_anomaly_total")
    success = Symbol("tracking.success_rate_possible_pct")
    sort!(result, [anomaly, success]; rev=[true, false])
    return result
end

function _worst_realizations(problem_runs::DataFrame)::DataFrame
    nrow(problem_runs) == 0 && return DataFrame()
    anomaly = Symbol("association_health.identity_anomaly_total")
    success = Symbol("tracking.success_rate_possible_pct")
    error = Symbol("tracking.converged_mean_error_m")
    rows = DataFrame[]
    for case_group in groupby(problem_runs, :case)
        ordered = sort(
            case_group,
            [anomaly, success, error];
            rev=[true, false, true]
        )
        push!(rows, DataFrame(ordered[1:1, :]))
    end
    return isempty(rows) ? DataFrame() : vcat(rows...; cols=:union)
end

function _difficulty_correlations(wide::DataFrame)::DataFrame
    nrow(wide) == 0 && return DataFrame()
    difficulties = (
        "geometry.max_simultaneously_visible_targets",
        "geometry.minimum_angular_separation_deg",
        "cross_m2m.iod_groups_mixed_target",
        "track_assoc.consensus_group_mixed_target"
    )
    outcomes = (
        "meas_assoc.recall_pct", "track_assoc.tt_recall_pct",
        "track_lifecycle.id_switch_total", "track_lifecycle.fragment_excess_total",
        "tracking.success_rate_possible_pct", "tracking.converged_mean_error_m"
    )
    rows = NamedTuple[]
    for case_group in groupby(wide, :case), difficulty in difficulties, outcome in outcomes
        dcol, ocol = Symbol(difficulty), Symbol(outcome)
        dcol in propertynames(case_group) && ocol in propertynames(case_group) || continue
        pairs = [(Float64(d), Float64(o)) for (d, o) in zip(case_group[!, dcol], case_group[!, ocol])
                 if isfinite(Float64(d)) && isfinite(Float64(o))]
        length(pairs) >= 3 || continue
        push!(rows, (
            case=String(case_group.case[1]), difficulty_metric=difficulty,
            outcome_metric=outcome, sample_count=length(pairs),
            pearson_correlation=cor(first.(pairs), last.(pairs))
        ))
    end
    return DataFrame(rows)
end

function _initialization_vs_convergence(initializations::DataFrame, windows::DataFrame)::DataFrame
    (nrow(initializations) == 0 || nrow(windows) == 0) && return DataFrame()
    rows = NamedTuple[]
    windows_by_key = Dict{Tuple{Int, String, Int, Int}, Vector{DataFrameRow}}()
    for window in eachrow(windows)
        key = (Int(window.realization), String(window.case), Int(window.observer), Int(window.target))
        push!(get!(windows_by_key, key, DataFrameRow[]), window)
    end
    for track in eachrow(initializations)
        target_id = Int(track.first_target_id)
        target_id > 0 || continue
        initialization_t = Float64(track.filter_initialized_t_s)
        isfinite(initialization_t) || continue
        key = (Int(track.realization), String(track.case), Int(track.observer), target_id)
        candidate_windows = get(windows_by_key, key, DataFrameRow[])
        matching_window = findfirst(candidate_windows) do window
            Float64(window.window_start_t_s) <= initialization_t <= Float64(window.window_end_t_s)
        end
        matching_window === nothing && continue
        window = candidate_windows[matching_window]
        push!(rows, (
            realization=Int(track.realization), case=String(track.case),
            scenario_seed=Int(track.scenario_seed), observer=Int(track.observer),
            slot=Int(track.slot), target=target_id, window_id=Int(window.window_id),
            initialization_t_s=initialization_t,
            initialization_latency_s=Float64(track.initialization_latency_s),
            initialization_position_error_m=Float64(track.initialization_position_error_m),
            converged_mean_error_m=Float64(window.converged_mean_error_m),
            converged_rmse_error_m=Float64(window.converged_rmse_error_m),
            iod_group_same_target=Int(track.iod_group_same_target),
            identity_class=String(track.identity_class),
            identity_switch_count=Int(track.id_switch_count),
            fragment_excess=Int(window.fragment_excess)
        ))
    end
    return DataFrame(rows)
end

function _write_aggregates(status::DataFrame)
    long = _read_all_metrics(status)
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_metrics_long.csv"), long)
    if nrow(long) == 0
        return nothing
    end
    long.metric_key = string.(long.section, ".", long.metric)
    wide = unstack(
        select(long, :realization, :case, :scenario_seed, :sensor_seed, :observer_od_seed, :metric_key, :value),
        [:realization, :case, :scenario_seed, :sensor_seed, :observer_od_seed],
        :metric_key,
        :value
    )
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_metrics_wide.csv"), wide)
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_summary.csv"), _metric_summary(long))
    problem_runs = _problem_characterization(wide)
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_problem_runs.csv"), problem_runs)
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_worst_runs.csv"), _worst_realizations(problem_runs))
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_difficulty_correlations.csv"), _difficulty_correlations(wide))
    tracking_windows = _combine_case_table(status, "tracking_window_table.csv")
    track_initializations = _combine_case_table(status, "track_initialization_table.csv")
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_tracking_windows.csv"), tracking_windows)
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_track_initializations.csv"), track_initializations)
    CSV.write(
        joinpath(MC_OUTPUT_ROOT, "mc_initialization_vs_convergence.csv"),
        _initialization_vs_convergence(track_initializations, tracking_windows)
    )
    CSV.write(
        joinpath(MC_OUTPUT_ROOT, "mc_target_populations.csv"),
        _combine_case_table(status, "target_population_table.csv"; one_case_per_realization=true)
    )
    CSV.write(
        joinpath(MC_OUTPUT_ROOT, "mc_iod_diagnostics.csv"),
        _combine_case_table(status, "iod_diagnostics_table.csv")
    )
    iod_pairwise = _combine_case_table(
        status,
        "iod_pairwise_consistency_table.csv"
    )
    if MC_SAVE_IOD_PAIRWISE || nrow(iod_pairwise) > 0
        CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_iod_pairwise.csv"), iod_pairwise)
    end

    fixed_names = Set((
        "m2t_gate_d2", "m2t_ratio_min", "t2t_iod_gate_d2",
        "t2t_filter_gate_d2", "t2t_ratio_min",
        "iod_position_rms_gate_m", "iod_validation_gate_d2"
    ))
    fixed = long[(long.section .== "run_config") .& in.(long.metric, Ref(fixed_names)), :]
    fixed_summary = combine(groupby(fixed, [:case, :metric]),
        :value => minimum => :minimum,
        :value => maximum => :maximum,
        :value => (values -> length(unique(values))) => :unique_value_count
    )
    any(fixed_summary.unique_value_count .!= 1) && error("Association thresholds changed across Monte Carlo runs")
    CSV.write(joinpath(MC_OUTPUT_ROOT, "mc_fixed_association_parameters.csv"), fixed_summary)
    return nothing
end

function _print_campaign_summary(status::DataFrame)
    completed = count(value -> value in ("ok", "resumed"), status.status)
    failed = count(==("failed"), status.status)
    println()
    println("Monte Carlo summary: completed=$(completed), failed=$(failed)")

    long = _read_all_metrics(status)
    nrow(long) == 0 && return nothing
    summary = _metric_summary(long)
    selected_metrics = (
        ("meas_assoc", "commit_accuracy_pct", "M2T precision [%]"),
        ("meas_assoc", "recall_pct", "M2T recall [%]"),
        ("track_assoc", "tt_accuracy_pct_known_only", "T2T precision [%]"),
        ("track_assoc", "tt_recall_pct", "T2T recall [%]"),
        ("track_lifecycle", "initialization_latency_median_s", "Initialization latency [s]"),
        ("track_lifecycle", "initialization_coverage_pct", "Initialization coverage [%]"),
        ("track_lifecycle", "initialization_success_pct", "Correct initialization success [%]"),
        ("track_lifecycle", "initialized_unique_targets", "Initialized trackable targets [count/run]"),
        ("track_lifecycle", "correctly_initialized_unique_targets", "Correctly initialized trackable targets [count/run]"),
        ("track_lifecycle", "never_initialized_unique_targets", "Never initialized trackable targets [count/run]"),
        ("track_lifecycle", "never_correctly_initialized_unique_targets", "Never correctly initialized targets [count/run]"),
        ("track_lifecycle", "id_switch_total", "Identity switches [count/run]"),
        ("track_lifecycle", "fragment_excess_total", "Track fragments [excess/run]"),
        ("track_lifecycle", "fragmented_window_rate_pct", "Fragmented tracking windows [%]"),
        ("track_lifecycle", "initialization_position_error_median_m", "Initialization error [m]"),
        ("tracking", "converged_median_error_m", "Converged error [m]"),
        ("track_lifecycle", "good_track_count", "Good tracks [count/run]"),
        ("track_lifecycle", "good_track_filter_duration_mean_s", "Good-track mean duration [s]"),
        ("track_lifecycle", "good_track_filter_duration_median_s", "Good-track median duration [s]"),
        ("track_lifecycle", "bad_track_count", "Bad tracks [count/run]"),
        ("track_lifecycle", "bad_track_filter_duration_mean_s", "Bad-track mean duration [s]"),
        ("track_lifecycle", "bad_track_filter_duration_median_s", "Bad-track median duration [s]"),
        ("cross_m2m", "iod_groups_mixed_target", "Wrong IOD groups [count/run]"),
        ("cross_m2m", "iod_position_cov_gate_rejected_same_target", "True IOD rejected by covariance [count/run]"),
        ("cross_m2m", "iod_position_cov_gate_rejected_mixed_target", "Wrong IOD rejected by covariance [count/run]"),
        ("cross_m2m", "iod_validation_rejected_same_target", "True IOD rejected by next-step score [count/run]"),
        ("cross_m2m", "iod_validation_rejected_mixed_target", "Wrong IOD rejected by next-step score [count/run]"),
        ("cross_m2m", "iod_validation_confirmed_mixed_target", "Wrong IOD promoted [count/run]"),
        ("cross_m2m", "iod_promoted_same_target_error_median_m", "Promoted true-IOD error [m]"),
        ("cross_m2m", "iod_promoted_mixed_target_error_median_m", "Promoted wrong-IOD error [m]"),
        ("track_assoc", "consensus_group_mixed_target", "Wrong consensus groups [count/run]"),
        ("tracking", "possible_windows", "Possible tracking windows [count/run]"),
        ("tracking", "tracked_windows", "Tracked windows [count/run]"),
        ("tracking", "successful_windows_under_1km", "Successful windows <1 km [count/run]"),
        ("tracking", "success_rate_possible_pct", "Tracking success [%]"),
        ("object_coverage", "jointly_detected_unique_targets", "Trackable unique targets [count/run]"),
        ("object_coverage", "successful_tracked_unique_targets", "Correctly tracked unique targets [count/run]"),
        ("object_coverage", "successful_unique_over_jointly_detected_pct", "Correctly tracked/trackable targets [%]"),
        ("geometry", "max_simultaneously_visible_targets", "Max simultaneously visible [targets]"),
        ("geometry", "minimum_angular_separation_deg", "Minimum angular separation [deg]")
    )

    for case in unique(String.(summary.case))
        println("  case=$(case)")
        for (section, metric, label) in selected_metrics
            rows = summary[
                (summary.case .== case) .&
                (summary.section .== section) .&
                (summary.metric .== metric),
                :
            ]
            nrow(rows) == 1 || continue
            row = rows[1, :]
            println(
                "    $(label): median=$(round(row.median; sigdigits=5)), " *
                "mean=$(round(row.mean; sigdigits=5)), " *
                "p05--p95=[$(round(row.p05; sigdigits=5)), $(round(row.p95; sigdigits=5))], " *
                "N=$(row.count)"
            )
        end
    end
    return nothing
end

function main()
    MC_RUNS > 0 || error("SPACEAGORA_MC_RUNS must be positive")
    cases = _selected_cases()
    mkpath(MC_OUTPUT_ROOT)
    status_path = joinpath(MC_OUTPUT_ROOT, "mc_run_status.csv")
    status = MC_RESUME && isfile(status_path) ?
        CSV.read(status_path, DataFrame; stringtype=String, pool=false) :
        DataFrame()

    println("Monte Carlo navigation campaign")
    println("  realizations: $(MC_FIRST_RUN):$(MC_FIRST_RUN + MC_RUNS - 1)")
    println("  cases: $(cases)")
    println("  observers: 16 (scenario definition)")
    println("  targets: $(get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", DEFAULT_TARGET_COUNT))")
    println("  mission time: $(get(ENV, "SPACEAGORA_MISSION_TIME_SEC", DEFAULT_MISSION_TIME_SEC)) s")
    println("  resume: $(MC_RESUME)")
    println("  save IOD event geometry: $(MC_SAVE_IOD_GEOMETRY)")
    println("  save IOD pairwise diagnostics: $(MC_SAVE_IOD_PAIRWISE)")
    println("  run IOD geometry summary: $(MC_IOD_GEOMETRY_ANALYSIS)")
    println("  output: $(MC_OUTPUT_ROOT)")
    flush(stdout)

    completed_cases = 0
    # Method-major ordering: complete every realization of one method before
    # moving to the next. This keeps partial campaign results easy to inspect.
    for case in cases, realization in MC_FIRST_RUN:(MC_FIRST_RUN + MC_RUNS - 1)
        row = _run_one(realization, case)
        status = _upsert_status!(status, row)
        CSV.write(status_path, status)
        completed_cases += 1
        if MC_AGGREGATE_EVERY > 0 && completed_cases % MC_AGGREGATE_EVERY == 0
            _write_aggregates(status)
            _print_campaign_summary(status)
        end
        row.status == "failed" && println(
            "  failed; inspect $(resolve_output_path(String(row.log_path), MC_OUTPUT_ROOT))"
        )
        flush(stdout)
    end
    _write_aggregates(status)
    if MC_AGGREGATE_EVERY <= 0 || completed_cases % MC_AGGREGATE_EVERY != 0
        _print_campaign_summary(status)
    end
    if MC_IOD_GEOMETRY_ANALYSIS
        geometry_case = get(
            ENV,
            "SPACEAGORA_MC_IOD_GEOMETRY_CASE",
            "proposed" in cases ? "proposed" : first(cases)
        )
        MonteCarloIODGeometryAnalysis.run_iod_geometry_analysis(
            MC_OUTPUT_ROOT;
            case=geometry_case,
            replay_missing=MC_IOD_GEOMETRY_REPLAY,
            max_realizations=MC_IOD_GEOMETRY_MAX_RUNS
        )
    end
    println("Campaign outputs written to $(MC_OUTPUT_ROOT)")
    flush(stdout)
end

main()
