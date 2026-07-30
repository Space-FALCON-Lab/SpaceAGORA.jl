using CSV
using DataFrames
using Statistics

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: REPO_ROOT, env_path, resolve_output_path, stored_output_path

const SIM_FILE = joinpath(@__DIR__, "..", "..", "navigation.jl")
const OUTPUT_ROOT = env_path(
    "SPACEAGORA_STRESS_OUTPUT",
    "sensitivity"
)
const N_RUNS = parse(Int, get(ENV, "SPACEAGORA_STRESS_RUNS", "100"))
const FIRST_RUN = parse(Int, get(ENV, "SPACEAGORA_STRESS_FIRST_RUN", "1"))
const BASE_SEED = parse(Int, get(ENV, "SPACEAGORA_STRESS_BASE_SEED", "202607180"))
const SEED_STRIDE = parse(Int, get(ENV, "SPACEAGORA_STRESS_SEED_STRIDE", "10"))
const RESUME = lowercase(strip(get(ENV, "SPACEAGORA_STRESS_RESUME", "true"))) in
    ("1", "true", "yes", "y", "on")
const STOP_ON_FAILURE = lowercase(strip(get(ENV, "SPACEAGORA_STRESS_STOP_ON_FAILURE", "true"))) in
    ("1", "true", "yes", "y", "on")
const AGGREGATE_EVERY = parse(Int, get(ENV, "SPACEAGORA_STRESS_AGGREGATE_EVERY", "100"))
const METRICS_SCHEMA_VERSION = 5
const DEFAULT_MISSION_TIME_SEC = "10000.0"
const DEFAULT_TARGET_COUNT = "300"
const DEFAULT_SIGMA_THETA_RAD = "1e-4"
const DEFAULT_HEAP_SIZE_HINT = "4G"
const DEFAULT_HEAVY_HEAP_SIZE_HINT = "4G"

Base.@kwdef struct StressSpec
    config::String
    sweep::String
    level::String
    stress_value::Float64
    stress_unit::String
    enable_bias::Bool = false
    bias_rad::Float64 = 0.0
    enable_misdetection::Bool = false
    misdetection_rate::Float64 = 0.0
    enable_false_alarm::Bool = false
    false_alarm_lambda::Float64 = 0.0
end

@inline _truthy(value::AbstractString)::Bool =
    lowercase(strip(value)) in ("1", "true", "yes", "y", "on")

const SAVE_IOD_EVENT_GEOMETRY = _truthy(get(
    ENV,
    "SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY",
    "false"
))
const SAVE_IOD_PAIRWISE_DIAGNOSTICS = _truthy(get(
    ENV,
    "SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS",
    "false"
))

function _selected_cases()::Vector{String}
    raw = get(ENV, "SPACEAGORA_STRESS_CASES", "proposed")
    cases = unique([strip(value) for value in split(raw, ',') if !isempty(strip(value))])
    valid = Set(("proposed", "centralized_oracle"))
    isempty(cases) && error("SPACEAGORA_STRESS_CASES must contain at least one case")
    invalid = [case for case in cases if !(case in valid)]
    isempty(invalid) || error("Invalid SPACEAGORA_STRESS_CASES entries: $(invalid)")
    # Preserve the requested order but always complete the proposed method first.
    sort!(cases; by=case -> case == "proposed" ? 0 : 1)
    return cases
end

function _selected_sweeps()::Set{String}
    raw = get(ENV, "SPACEAGORA_STRESS_SWEEPS", "bias,misdetection,false_alarm")
    sweeps = Set(strip(value) for value in split(raw, ',') if !isempty(strip(value)))
    valid = Set(("bias", "misdetection", "false_alarm"))
    isempty(sweeps) && error("SPACEAGORA_STRESS_SWEEPS must contain at least one sweep")
    invalid = setdiff(sweeps, valid)
    isempty(invalid) || error("Invalid SPACEAGORA_STRESS_SWEEPS entries: $(collect(invalid))")
    return sweeps
end

function _stress_specs()::Vector{StressSpec}
    selected = _selected_sweeps()
    specs = StressSpec[
        StressSpec(
            config="nominal", sweep="nominal", level="nominal",
            stress_value=0.0, stress_unit="none"
        )
    ]
    if "bias" in selected
        for (level, value) in zip(("light", "moderate", "stress"), (1e-5, 5e-5, 1e-4))
            push!(specs, StressSpec(
                config="bias_$(level)", sweep="bias", level=level,
                stress_value=value, stress_unit="rad",
                enable_bias=true, bias_rad=value
            ))
        end
    end
    if "misdetection" in selected
        for (level, value) in zip(("light", "moderate", "stress"), (0.05, 0.10, 0.20))
            push!(specs, StressSpec(
                config="misdetection_$(level)", sweep="misdetection", level=level,
                stress_value=value, stress_unit="probability",
                enable_misdetection=true, misdetection_rate=value
            ))
        end
    end
    if "false_alarm" in selected
        for (level, value) in zip(
            ("stress", "severe", "extreme"),
            (0.20, 0.40, 0.80)
        )
            push!(specs, StressSpec(
                config="false_alarm_$(level)", sweep="false_alarm", level=level,
                stress_value=value, stress_unit="mean_per_observer_epoch",
                enable_false_alarm=true, false_alarm_lambda=value
            ))
        end
    end
    return specs
end

@inline function _run_seeds(realization::Int)
    scenario = BASE_SEED + SEED_STRIDE * (realization - 1)
    return (
        scenario=scenario,
        sensor=scenario + 1,
        od=scenario + 2,
        bias=scenario + 3,
        misdetection=scenario + 4,
        false_alarm_count=scenario + 5,
        false_alarm_direction=scenario + 6
    )
end

@inline _run_dir(realization::Int) =
    joinpath(OUTPUT_ROOT, "run_$(lpad(realization, 4, '0'))")
@inline _config_dir(realization::Int, spec::StressSpec) =
    joinpath(_run_dir(realization), spec.config)
@inline _case_dir(realization::Int, spec::StressSpec, case::String) =
    joinpath(_config_dir(realization, spec), case)
@inline _metrics_path(realization::Int, spec::StressSpec, case::String) =
    joinpath(_case_dir(realization, spec, case), "association_quality_table.csv")

const REQUIRED_CASE_OUTPUTS = (
    "association_quality_table.csv",
)

function _outputs_complete(case_dir::String, log_path::String)::Bool
    all(name -> begin
        output_path = joinpath(case_dir, name)
        isfile(output_path) && filesize(output_path) > 0
    end, REQUIRED_CASE_OUTPUTS) || return false
    isfile(log_path) || return false
    log_text = read(log_path, String)
    # The aggregate table is written atomically before its checkpoint. Accept
    # that checkpoint as a legacy completion marker for runs produced before
    # the stress campaign switched to metrics-only output.
    return occursin("CSV_WRITE_DONE", log_text) ||
        occursin("stage=csv_association_quality", log_text)
end

function _metric_value(table::DataFrame, section::String, metric::String)
    rows = table[(table.section .== section) .& (table.metric .== metric), :]
    return nrow(rows) == 1 ? Float64(rows.value[1]) : NaN
end

function _completed(realization::Int, spec::StressSpec, case::String, seeds)::Bool
    path = _metrics_path(realization, spec, case)
    case_dir = _case_dir(realization, spec, case)
    log_path = joinpath(case_dir, "run.log")
    _outputs_complete(case_dir, log_path) || return false
    try
        table = CSV.read(path, DataFrame)
        checks = (
            _metric_value(table, "run_config", "scenario_seed") == seeds.scenario,
            _metric_value(table, "run_config", "sensor_seed") == seeds.sensor,
            _metric_value(table, "run_config", "observer_od_seed") == seeds.od,
            _metric_value(table, "run_config", "mc_metrics_schema_version") == METRICS_SCHEMA_VERSION,
            !spec.enable_misdetection ||
                _metric_value(table, "run_config", "tracking_windows_truth_visibility") == 1.0,
            _metric_value(table, "sensor_errors", "enable_measurement_bias") == (spec.enable_bias ? 1.0 : 0.0),
            _metric_value(table, "sensor_errors", "enable_missed_detections") == (spec.enable_misdetection ? 1.0 : 0.0),
            _metric_value(table, "sensor_errors", "enable_false_alarms") == (spec.enable_false_alarm ? 1.0 : 0.0),
            !SAVE_IOD_PAIRWISE_DIAGNOSTICS ||
                isfile(joinpath(case_dir, "iod_pairwise_consistency_table.csv")),
            isapprox(_metric_value(table, "sensor_errors", "configured_measurement_bias_rad"), spec.bias_rad),
            isapprox(_metric_value(table, "sensor_errors", "configured_missed_detection_rate"), spec.misdetection_rate),
            isapprox(_metric_value(table, "sensor_errors", "configured_false_alarm_rate"), spec.false_alarm_lambda)
        )
        return all(checks)
    catch
        return false
    end
end

function _child_environment(realization::Int, spec::StressSpec, case::String, seeds)
    child_env = copy(ENV)
    child_env["SPACEAGORA_COMPARISON_OUTPUT"] = _config_dir(realization, spec)
    child_env["SPACEAGORA_NAV_CASE"] = case
    child_env["SPACEAGORA_NAV_OUTPUT_LABEL"] = case
    child_env["SPACEAGORA_SCENARIO_SEED"] = string(seeds.scenario)
    child_env["SPACEAGORA_SENSOR_SEED"] = string(seeds.sensor)
    child_env["SPACEAGORA_OBSERVER_OD_SEED"] = string(seeds.od)
    child_env["SPACEAGORA_MEASUREMENT_BIAS_SEED"] = string(seeds.bias)
    child_env["SPACEAGORA_MISDETECTION_SEED"] = string(seeds.misdetection)
    child_env["SPACEAGORA_FALSE_ALARM_COUNT_SEED"] = string(seeds.false_alarm_count)
    child_env["SPACEAGORA_FALSE_ALARM_DIRECTION_SEED"] = string(seeds.false_alarm_direction)
    child_env["SPACEAGORA_MISSION_TIME_SEC"] =
        get(ENV, "SPACEAGORA_MISSION_TIME_SEC", DEFAULT_MISSION_TIME_SEC)
    child_env["SPACEAGORA_N_CLUSTER_TARGETS"] =
        get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", DEFAULT_TARGET_COUNT)
    child_env["SPACEAGORA_SIGMA_THETA_RAD"] =
        get(ENV, "SPACEAGORA_SIGMA_THETA_RAD", DEFAULT_SIGMA_THETA_RAD)
    child_env["SPACEAGORA_ENABLE_MEASUREMENT_BIAS"] = string(spec.enable_bias)
    child_env["SPACEAGORA_MEASUREMENT_BIAS_RAD"] = string(spec.bias_rad)
    child_env["SPACEAGORA_ENABLE_MISDETECTIONS"] = string(spec.enable_misdetection)
    child_env["SPACEAGORA_MISDETECTION_RATE"] = string(spec.misdetection_rate)
    child_env["SPACEAGORA_ENABLE_FALSE_ALARMS"] = string(spec.enable_false_alarm)
    child_env["SPACEAGORA_FALSE_ALARM_RATE"] = string(spec.false_alarm_lambda)
    child_env["SPACEAGORA_SAVE_SIMULATION_RESULTS"] = "false"
    child_env["SPACEAGORA_SAVE_TARGET_ESTIMATE_FIELDS"] = "false"
    child_env["SPACEAGORA_SAVE_BUNDLE"] = "0"
    child_env["SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES"] =
        get(ENV, "SPACEAGORA_STRESS_SAVE_DETAILED_TABLES", "false")
    child_env["SPACEAGORA_SAVE_AUXILIARY_METRIC_TABLES"] =
        get(ENV, "SPACEAGORA_STRESS_SAVE_AUXILIARY_METRIC_TABLES", "false")
    child_env["SPACEAGORA_ENABLE_NAV_TIMING"] =
        get(ENV, "SPACEAGORA_STRESS_ENABLE_NAV_TIMING", "false")
    child_env["SPACEAGORA_SAVE_IOD_EVENT_GEOMETRY"] =
        string(SAVE_IOD_EVENT_GEOMETRY)
    child_env["SPACEAGORA_SAVE_IOD_PAIRWISE_DIAGNOSTICS"] =
        string(SAVE_IOD_PAIRWISE_DIAGNOSTICS)
    child_env["SPACEAGORA_LOG_NAV_EVENTS"] = "false"
    child_env["SPACEAGORA_LOG_SIM_PROGRESS"] =
        get(ENV, "SPACEAGORA_STRESS_LOG_SIM_PROGRESS", "true")
    child_env["SPACEAGORA_SIM_PROGRESS_INTERVAL_SEC"] =
        get(ENV, "SPACEAGORA_STRESS_PROGRESS_INTERVAL_SEC", "1000.0")
    child_env["SPACEAGORA_FORCE_INCREMENTAL_GC"] =
        get(ENV, "SPACEAGORA_STRESS_FORCE_INCREMENTAL_GC", "true")
    return child_env
end

function _heap_size_hint(spec::StressSpec)::String
    # Heavy cases have a dedicated override so a global 4G setting cannot
    # disable the earlier garbage collection used by this safeguard.
    memory_heavy = spec.stress_value >= 0.10 &&
        spec.sweep in ("misdetection", "false_alarm")
    return memory_heavy ?
        strip(get(ENV, "SPACEAGORA_STRESS_HEAVY_HEAP_SIZE_HINT", DEFAULT_HEAVY_HEAP_SIZE_HINT)) :
        strip(get(ENV, "SPACEAGORA_JULIA_HEAP_SIZE_HINT", DEFAULT_HEAP_SIZE_HINT))
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
    error_index !== nothing &&
        return strip(join(lines[error_index:min(error_index + 2, length(lines))], " | "))
    tail = [strip(line) for line in lines[max(1, length(lines) - 4):end] if !isempty(strip(line))]
    return isempty(tail) ? "run log contains only blank lines" : join(tail, " | ")
end

function _status_row(realization::Int, spec::StressSpec, case::String, status::String,
                     exit_code::Int, elapsed_sec::Float64, seeds, case_dir::String,
                     log_path::String)
    return (
        realization=realization, case=case, config=spec.config, sweep=spec.sweep,
        level=spec.level, stress_value=spec.stress_value, stress_unit=spec.stress_unit,
        enable_bias=spec.enable_bias ? 1 : 0, bias_rad=spec.bias_rad,
        enable_misdetection=spec.enable_misdetection ? 1 : 0,
        misdetection_rate=spec.misdetection_rate,
        enable_false_alarm=spec.enable_false_alarm ? 1 : 0,
        false_alarm_lambda=spec.false_alarm_lambda,
        status=status, exit_code=exit_code, elapsed_sec=elapsed_sec,
        scenario_seed=seeds.scenario, sensor_seed=seeds.sensor,
        observer_od_seed=seeds.od, measurement_bias_seed=seeds.bias,
        misdetection_seed=seeds.misdetection,
        false_alarm_count_seed=seeds.false_alarm_count,
        false_alarm_direction_seed=seeds.false_alarm_direction,
        results_dir=stored_output_path(case_dir, OUTPUT_ROOT),
        log_path=stored_output_path(log_path, OUTPUT_ROOT)
    )
end

function _run_one(realization::Int, spec::StressSpec, case::String)
    seeds = _run_seeds(realization)
    case_dir = _case_dir(realization, spec, case)
    mkpath(case_dir)
    log_path = joinpath(case_dir, "run.log")
    if RESUME && _completed(realization, spec, case, seeds)
        println("  resumed run=$(realization), config=$(spec.config), case=$(case)")
        return _status_row(realization, spec, case, "resumed", 0, 0.0, seeds, case_dir, log_path)
    end

    # A previously interrupted child may have left an older metrics file in
    # the case directory.  Remove the completion marker before retrying so
    # that success can only be attributed to the current child process.
    metrics_path = _metrics_path(realization, spec, case)
    isfile(metrics_path) && rm(metrics_path; force=true)

    start_time = time()
    julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
    heap_hint = _heap_size_hint(spec)
    child_env = _child_environment(realization, spec, case, seeds)
    cmd = isempty(heap_hint) ?
        setenv(`$julia_exe --project=$REPO_ROOT $SIM_FILE`, child_env) :
        setenv(`$julia_exe --heap-size-hint=$heap_hint --project=$REPO_ROOT $SIM_FILE`, child_env)
    println("  run=$(realization), config=$(spec.config), case=$(case), " *
            "seed=$(seeds.scenario), heap_hint=$(isempty(heap_hint) ? "none" : heap_hint)")
    flush(stdout)
    exit_code = 1
    open(log_path, "w") do log_io
        process = run(pipeline(ignorestatus(cmd); stdout=log_io, stderr=log_io); wait=false)
        heartbeat = time()
        while process_running(process)
            sleep(1.0)
            if time() - heartbeat >= 30.0
                println("    still running ($(round((time() - start_time) / 60; digits=1)) min)")
                flush(stdout)
                heartbeat = time()
            end
        end
        wait(process)
        exit_code = try Int(process.exitcode) catch; success(process) ? 0 : 1 end
    end
    outputs_complete = _outputs_complete(case_dir, log_path)
    ok = exit_code == 0 && outputs_complete
    elapsed = time() - start_time
    if ok
        println("    completed in $(round(elapsed; digits=1)) s")
    else
        reason = exit_code != 0 ? "child exit code $(exit_code)" :
            "child ended before the required stress metrics were committed"
        println("    FAILED ($(reason)): $(_failure_excerpt(log_path))")
    end
    flush(stdout)
    return _status_row(
        realization, spec, case, ok ? "ok" : "failed", exit_code, elapsed,
        seeds, case_dir, log_path
    )
end

function _upsert_status!(status::DataFrame, row)::DataFrame
    ncol(status) == 0 && return DataFrame([row])
    matching_rows = findall(
        (status.realization .== row.realization) .&
        (status.case .== row.case) .&
        (status.config .== row.config)
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
        String(row.status) in ("ok", "resumed") || continue
        results_dir = resolve_output_path(String(row.results_dir), OUTPUT_ROOT)
        path = joinpath(results_dir, "association_quality_table.csv")
        isfile(path) || continue
        table = CSV.read(path, DataFrame)
        insertcols!(table, 1,
            :realization => fill(Int(row.realization), nrow(table)),
            :case => fill(String(row.case), nrow(table)),
            :config => fill(String(row.config), nrow(table)),
            :sweep => fill(String(row.sweep), nrow(table)),
            :level => fill(String(row.level), nrow(table)),
            :stress_value => fill(Float64(row.stress_value), nrow(table)),
            :stress_unit => fill(String(row.stress_unit), nrow(table)),
            :scenario_seed => fill(Int(row.scenario_seed), nrow(table)),
            :sensor_seed => fill(Int(row.sensor_seed), nrow(table)),
            :observer_od_seed => fill(Int(row.observer_od_seed), nrow(table))
        )
        push!(tables, table)
    end
    return isempty(tables) ? DataFrame() : vcat(tables...; cols=:union)
end

function _combine_table(
    status::DataFrame,
    filename::String;
    one_per_realization::Bool=false
)::DataFrame
    tables = DataFrame[]
    seen_realizations = Set{Int}()
    for row in eachrow(status)
        String(row.status) in ("ok", "resumed") || continue
        realization = Int(row.realization)
        one_per_realization && realization in seen_realizations && continue
        results_dir = resolve_output_path(String(row.results_dir), OUTPUT_ROOT)
        path = joinpath(results_dir, filename)
        isfile(path) || continue
        table = CSV.read(path, DataFrame)
        insertcols!(table, 1,
            :realization => fill(realization, nrow(table)),
            :case => fill(String(row.case), nrow(table)),
            :config => fill(String(row.config), nrow(table)),
            :sweep => fill(String(row.sweep), nrow(table)),
            :level => fill(String(row.level), nrow(table)),
            :stress_value => fill(Float64(row.stress_value), nrow(table)),
            :scenario_seed => fill(Int(row.scenario_seed), nrow(table))
        )
        push!(tables, table)
        one_per_realization && push!(seen_realizations, realization)
    end
    return isempty(tables) ? DataFrame() : vcat(tables...; cols=:union)
end

function _metric_summary(long::DataFrame)::DataFrame
    rows = NamedTuple[]
    nrow(long) == 0 && return DataFrame(rows)
    keys = [:case, :config, :sweep, :level, :stress_value, :stress_unit, :section, :metric]
    for group in groupby(long, keys)
        values = Float64[Float64(value) for value in group.value if isfinite(Float64(value))]
        isempty(values) && continue
        push!(rows, (
            case=String(group.case[1]), config=String(group.config[1]),
            sweep=String(group.sweep[1]), level=String(group.level[1]),
            stress_value=Float64(group.stress_value[1]), stress_unit=String(group.stress_unit[1]),
            section=String(group.section[1]), metric=String(group.metric[1]),
            count=length(values), mean=mean(values),
            std=length(values) > 1 ? std(values) : NaN,
            minimum=minimum(values), p05=quantile(values, 0.05),
            median=median(values), p95=quantile(values, 0.95), maximum=maximum(values)
        ))
    end
    return DataFrame(rows)
end

function _paired_deltas(long::DataFrame)::DataFrame
    nrow(long) == 0 && return DataFrame()
    nominal = long[long.config .== "nominal", :]
    stressed = long[long.config .!= "nominal", :]
    nrow(nominal) == 0 && return DataFrame()
    nominal = select(
        nominal, :realization, :case, :section, :metric,
        :value => :nominal_value
    )
    paired = innerjoin(stressed, nominal; on=[:realization, :case, :section, :metric])
    paired.delta = paired.value .- paired.nominal_value
    return paired
end

function _problem_characterization(wide::DataFrame)::DataFrame
    nrow(wide) == 0 && return DataFrame()
    wanted = [
        "sensor_errors.realized_bias_component_rms_rad",
        "sensor_errors.realized_missed_rate_pct",
        "sensor_errors.realized_false_alarm_mean_per_observer_epoch",
        "association_health.identity_anomaly_total",
        "meas_assoc.commit_accuracy_pct", "meas_assoc.recall_pct",
        "track_assoc.tt_accuracy_pct_known_only", "track_assoc.tt_recall_pct",
        "cross_m2m.iod_groups_mixed_target",
        "cross_m2m.iod_validation_confirmed_mixed_target",
        "track_assoc.consensus_group_mixed_target",
        "false_alarm.fake_tracks_initialized",
        "track_lifecycle.id_switch_total", "track_lifecycle.fragment_excess_total",
        "track_lifecycle.initialization_success_pct",
        "track_lifecycle.never_correctly_initialized_unique_targets",
        "track_lifecycle.initialization_position_error_mean_m",
        "tracking.tracking_coverage_pct", "tracking.success_rate_possible_pct",
        "tracking.mean_error_tracked_windows_m", "tracking.converged_mean_error_m",
        "object_coverage.successful_unique_over_jointly_detected_pct",
        "geometry.max_simultaneously_visible_targets",
        "geometry.minimum_angular_separation_deg"
    ]
    columns = [
        :realization, :case, :config, :sweep, :level, :stress_value,
        :scenario_seed, :sensor_seed, :observer_od_seed
    ]
    append!(columns, [Symbol(name) for name in wanted if Symbol(name) in propertynames(wide)])
    result = select(wide, columns)
    anomaly = Symbol("association_health.identity_anomaly_total")
    success = Symbol("tracking.success_rate_possible_pct")
    sort!(result, [:case, :config, anomaly, success]; rev=[false, false, true, false])
    return result
end

function _worst_realizations(problem_runs::DataFrame)::DataFrame
    nrow(problem_runs) == 0 && return DataFrame()
    anomaly = Symbol("association_health.identity_anomaly_total")
    success = Symbol("tracking.success_rate_possible_pct")
    error = Symbol("tracking.converged_mean_error_m")
    rows = DataFrame[]
    for group in groupby(problem_runs, [:case, :config])
        ordered = sort(group, [anomaly, success, error]; rev=[true, false, true])
        push!(rows, DataFrame(ordered[1:1, :]))
    end
    return isempty(rows) ? DataFrame() : vcat(rows...; cols=:union)
end

function _difficulty_correlations(wide::DataFrame)::DataFrame
    nrow(wide) == 0 && return DataFrame()
    difficulties = (
        "sensor_errors.realized_bias_component_rms_rad",
        "sensor_errors.realized_missed_rate_pct",
        "sensor_errors.realized_false_alarm_mean_per_observer_epoch",
        "geometry.max_simultaneously_visible_targets",
        "geometry.minimum_angular_separation_deg",
        "cross_m2m.iod_groups_mixed_target",
        "track_assoc.consensus_group_mixed_target"
    )
    outcomes = (
        "meas_assoc.recall_pct", "track_assoc.tt_recall_pct",
        "track_lifecycle.initialization_success_pct",
        "track_lifecycle.id_switch_total", "track_lifecycle.fragment_excess_total",
        "tracking.tracking_coverage_pct", "tracking.success_rate_possible_pct",
        "tracking.converged_mean_error_m"
    )
    rows = NamedTuple[]
    for group in groupby(wide, [:case, :config]), difficulty in difficulties, outcome in outcomes
        dcol, ocol = Symbol(difficulty), Symbol(outcome)
        dcol in propertynames(group) && ocol in propertynames(group) || continue
        pairs = [
            (Float64(d), Float64(o))
            for (d, o) in zip(group[!, dcol], group[!, ocol])
            if isfinite(Float64(d)) && isfinite(Float64(o))
        ]
        length(pairs) >= 3 || continue
        push!(rows, (
            case=String(group.case[1]), config=String(group.config[1]),
            sweep=String(group.sweep[1]), level=String(group.level[1]),
            stress_value=Float64(group.stress_value[1]),
            difficulty_metric=difficulty, outcome_metric=outcome,
            sample_count=length(pairs),
            pearson_correlation=cor(first.(pairs), last.(pairs))
        ))
    end
    return DataFrame(rows)
end

function _initialization_vs_convergence(
    initializations::DataFrame,
    windows::DataFrame
)::DataFrame
    (nrow(initializations) == 0 || nrow(windows) == 0) && return DataFrame()
    rows = NamedTuple[]
    key_type = Tuple{Int, String, String, Int, Int}
    windows_by_key = Dict{key_type, Vector{DataFrameRow}}()
    for window in eachrow(windows)
        key = (
            Int(window.realization), String(window.case), String(window.config),
            Int(window.observer), Int(window.target)
        )
        push!(get!(windows_by_key, key, DataFrameRow[]), window)
    end
    for track in eachrow(initializations)
        target_id = Int(track.first_target_id)
        target_id > 0 || continue
        initialization_t = Float64(track.filter_initialized_t_s)
        isfinite(initialization_t) || continue
        key = (
            Int(track.realization), String(track.case), String(track.config),
            Int(track.observer), target_id
        )
        candidate_windows = get(windows_by_key, key, DataFrameRow[])
        matching_window = findfirst(candidate_windows) do window
            Float64(window.window_start_t_s) <= initialization_t <=
                Float64(window.window_end_t_s)
        end
        matching_window === nothing && continue
        window = candidate_windows[matching_window]
        push!(rows, (
            realization=Int(track.realization), case=String(track.case),
            config=String(track.config), sweep=String(track.sweep),
            level=String(track.level), stress_value=Float64(track.stress_value),
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
    CSV.write(joinpath(OUTPUT_ROOT, "stress_metrics_long.csv"), long)
    nrow(long) == 0 && return nothing
    long.metric_key = string.(long.section, ".", long.metric)
    identifiers = [
        :realization, :case, :config, :sweep, :level, :stress_value,
        :stress_unit, :scenario_seed, :sensor_seed, :observer_od_seed
    ]
    wide = unstack(select(long, identifiers..., :metric_key, :value), identifiers, :metric_key, :value)
    CSV.write(joinpath(OUTPUT_ROOT, "stress_metrics_wide.csv"), wide)
    CSV.write(joinpath(OUTPUT_ROOT, "stress_summary.csv"), _metric_summary(long))
    problem_runs = _problem_characterization(wide)
    CSV.write(joinpath(OUTPUT_ROOT, "stress_problem_runs.csv"), problem_runs)
    CSV.write(joinpath(OUTPUT_ROOT, "stress_worst_runs.csv"), _worst_realizations(problem_runs))
    CSV.write(
        joinpath(OUTPUT_ROOT, "stress_difficulty_correlations.csv"),
        _difficulty_correlations(wide)
    )

    paired = _paired_deltas(long)
    CSV.write(joinpath(OUTPUT_ROOT, "stress_paired_deltas.csv"), paired)
    if nrow(paired) > 0
        paired_for_summary = select(paired, Not(:value))
        rename!(paired_for_summary, :delta => :value)
        CSV.write(
            joinpath(OUTPUT_ROOT, "stress_paired_delta_summary.csv"),
            _metric_summary(paired_for_summary)
        )
    end

    tracking_windows = _combine_table(status, "tracking_window_table.csv")
    track_initializations = _combine_table(status, "track_initialization_table.csv")
    CSV.write(joinpath(OUTPUT_ROOT, "stress_tracking_windows.csv"), tracking_windows)
    CSV.write(joinpath(OUTPUT_ROOT, "stress_track_initializations.csv"), track_initializations)
    CSV.write(
        joinpath(OUTPUT_ROOT, "stress_initialization_vs_convergence.csv"),
        _initialization_vs_convergence(track_initializations, tracking_windows)
    )
    CSV.write(
        joinpath(OUTPUT_ROOT, "stress_iod_diagnostics.csv"),
        _combine_table(status, "iod_diagnostics_table.csv")
    )
    iod_pairwise = _combine_table(
        status,
        "iod_pairwise_consistency_table.csv"
    )
    if SAVE_IOD_PAIRWISE_DIAGNOSTICS || nrow(iod_pairwise) > 0
        CSV.write(
            joinpath(OUTPUT_ROOT, "stress_iod_pairwise.csv"),
            iod_pairwise
        )
    end
    CSV.write(
        joinpath(OUTPUT_ROOT, "stress_target_populations.csv"),
        _combine_table(status, "target_population_table.csv"; one_per_realization=true)
    )

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
    any(fixed_summary.unique_value_count .!= 1) &&
        error("Association thresholds changed across stress configurations")
    CSV.write(joinpath(OUTPUT_ROOT, "stress_fixed_association_parameters.csv"), fixed_summary)
    return nothing
end

const COMMON_METRICS = (
    ("meas_assoc", "commit_accuracy_pct", "M2T precision [%]"),
    ("meas_assoc", "recall_pct", "M2T recall [%]"),
    ("track_assoc", "tt_accuracy_pct_known_only", "T2T precision [%]"),
    ("track_assoc", "tt_recall_pct", "T2T recall [%]"),
    ("track_lifecycle", "initialization_success_pct", "Initialization success [%]"),
    ("track_lifecycle", "initialization_latency_median_s", "Initialization latency [s]"),
    ("tracking", "success_rate_possible_pct", "Tracking success [%]"),
    ("tracking", "converged_median_error_m", "Converged error [m]"),
    ("track_lifecycle", "fragment_excess_total", "Track fragments [count/run]"),
    ("track_lifecycle", "id_switch_total", "Identity switches [count/run]")
)

const SWEEP_METRICS = Dict(
    "bias" => (
        ("sensor_errors", "realized_bias_component_rms_rad", "Realized bias RMS [rad]"),
        ("track_lifecycle", "initialization_position_error_median_m", "Initialization error [m]"),
        ("cross_m2m", "iod_validation_rejected_same_target", "True IOD rejected at validation [count/run]"),
        ("track_assoc", "consensus_group_mixed_target", "Wrong consensus groups [count/run]")
    ),
    "misdetection" => (
        ("sensor_errors", "realized_missed_rate_pct", "Realized missed detections [%]"),
        ("track_lifecycle", "never_correctly_initialized_unique_targets", "Never initialized correctly [count/run]"),
        ("cross_m2m", "iod_validation_no_measure", "IOD validations without measurement [count/run]"),
        ("track_lifecycle", "fragmented_window_rate_pct", "Fragmented windows [%]")
    ),
    "false_alarm" => (
        ("sensor_errors", "realized_false_alarm_mean_per_observer_epoch", "False alarms / observer-epoch"),
        ("sensor_errors", "realized_false_alarm_multiple_epoch_pct", "Observer-epochs with >=2 false alarms [%]"),
        ("false_alarm", "committed_to_m2t", "False alarms committed to M2T [count/run]"),
        ("false_alarm", "seed_tracks_created", "False-alarm seed tracks [count/run]"),
        ("false_alarm", "iod_groups_false_contaminated", "False-contaminated IOD groups [count/run]"),
        ("false_alarm", "iod_groups_false_only", "False-only IOD groups [count/run]"),
        ("false_alarm", "iod_validation_rejected_false_contaminated", "Contaminated IODs rejected [count/run]"),
        ("false_alarm", "iod_initialized_false_contaminated", "Contaminated IODs promoted [count/run]"),
        ("false_alarm", "fake_tracks_initialized", "Fake initialized tracks [count/run]"),
        ("false_alarm", "fake_track_duration_median_s", "Fake-track duration [s]")
    )
)

function _print_summary(status::DataFrame)
    completed = count(value -> String(value) in ("ok", "resumed"), status.status)
    failed = count(value -> String(value) == "failed", status.status)
    println("Stress Monte Carlo summary: completed=$(completed), failed=$(failed)")
    long = _read_all_metrics(status)
    nrow(long) == 0 && return nothing
    summary = _metric_summary(long)
    configs = unique(select(summary, :case, :config, :sweep, :level, :stress_value))
    sort!(configs, [:case, :sweep, :stress_value])
    for config_row in eachrow(configs)
        case, config, sweep = String(config_row.case), String(config_row.config), String(config_row.sweep)
        println("  case=$(case), config=$(config), value=$(config_row.stress_value)")
        metrics = sweep == "nominal" ? COMMON_METRICS : (COMMON_METRICS..., SWEEP_METRICS[sweep]...)
        for (section, metric, label) in metrics
            rows = summary[
                (summary.case .== case) .& (summary.config .== config) .&
                (summary.section .== section) .& (summary.metric .== metric), :
            ]
            nrow(rows) == 1 || continue
            row = rows[1, :]
            println(
                "    $(label): median=$(round(row.median; sigdigits=5)), " *
                "p05--p95=[$(round(row.p05; sigdigits=5)), $(round(row.p95; sigdigits=5))], " *
                "N=$(row.count)"
            )
        end
    end
    return nothing
end

function analyze_existing()::Nothing
    status_path = joinpath(OUTPUT_ROOT, "stress_run_status.csv")
    isfile(status_path) || error("Missing sensitivity Monte Carlo status: $(status_path)")
    status = CSV.read(status_path, DataFrame)
    active_configs = Set(spec.config for spec in _stress_specs())
    saved_configs = Set(String.(status.config))
    ignored_configs = sort!(collect(setdiff(saved_configs, active_configs)))
    if !isempty(ignored_configs)
        println(
            "Ignoring sensitivity configurations not used by the current " *
            "campaign: $(join(ignored_configs, ", "))"
        )
    end
    status = status[
        [String(config) in active_configs for config in status.config],
        :
    ]
    nrow(status) > 0 || error(
        "No saved runs match the currently selected sensitivity configurations."
    )
    _write_aggregates(status)
    _print_summary(status)
    println("Sensitivity Monte Carlo analysis written to $(OUTPUT_ROOT)")
    return nothing
end

function main()
    N_RUNS > 0 || error("SPACEAGORA_STRESS_RUNS must be positive")
    cases = _selected_cases()
    specs = _stress_specs()
    mkpath(OUTPUT_ROOT)
    status_path = joinpath(OUTPUT_ROOT, "stress_run_status.csv")
    status = RESUME && isfile(status_path) ?
        CSV.read(status_path, DataFrame; stringtype=String, pool=false) :
        DataFrame()
    total_jobs = length(cases) * length(specs) * N_RUNS
    println("Navigation stress Monte Carlo campaign")
    println("  realizations: $(FIRST_RUN):$(FIRST_RUN + N_RUNS - 1)")
    println("  cases: $(join(cases, ", ")) (proposed is always completed first)")
    println("  configurations: $(join(getfield.(specs, :config), ", "))")
    println("  total jobs: $(total_jobs)")
    println("  false alarms: Poisson count per observer-epoch")
    println("  common random numbers: enabled across stress levels")
    println("  save IOD event geometry: $(SAVE_IOD_EVENT_GEOMETRY)")
    println("  save IOD pairwise diagnostics: $(SAVE_IOD_PAIRWISE_DIAGNOSTICS)")
    println("  output: $(OUTPUT_ROOT)")
    flush(stdout)

    completed_jobs = 0
    # Method-major order makes the proposed campaign self-contained even if
    # the optional centralized comparison is interrupted or omitted.
    for case in cases, spec in specs, realization in FIRST_RUN:(FIRST_RUN + N_RUNS - 1)
        row = _run_one(realization, spec, case)
        status = _upsert_status!(status, row)
        CSV.write(status_path, status)
        completed_jobs += 1
        if String(row.status) == "failed" && STOP_ON_FAILURE
            error(
                "Stress campaign stopped after failed job: " *
                "run=$(row.realization), config=$(row.config), case=$(row.case). " *
                "Completed outputs remain resumable."
            )
        end
        if AGGREGATE_EVERY > 0 && completed_jobs % AGGREGATE_EVERY == 0
            _write_aggregates(status)
            # Aggregation temporarily materializes all completed compact CSV
            # tables. Release those parent-process buffers before spawning the
            # next independent simulation child.
            GC.gc(true)
            println("  aggregate checkpoint: $(completed_jobs)/$(total_jobs) jobs")
            flush(stdout)
        end
    end
    _write_aggregates(status)
    _print_summary(status)
    println("Campaign outputs written to $(OUTPUT_ROOT)")
    return nothing
end

if _truthy(get(ENV, "SPACEAGORA_STRESS_ANALYZE_ONLY", "false"))
    Base.invokelatest(analyze_existing)
else
    main()
end
