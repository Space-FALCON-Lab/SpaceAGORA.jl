using Dates
using CSV
using DataFrames
using Random
using Printf

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: REPO_ROOT, env_override_path, resolve_output_path
using .NavigationPaths: stored_output_path, transient_navigation_path

include(joinpath(@__DIR__, "compact_summary.jl"))

const SIM_FILE = joinpath(@__DIR__, "..", "..", "navigation.jl")
const OUTPUT_ROOT = env_override_path(
    "SPACEAGORA_SENSITIVITY_OUTPUT",
    transient_navigation_path("runs", "sensitivity_sweep")
)

const DEFAULT_SCENARIO_SEED = "1005639417"
const DEFAULT_START_UTC = "2026-06-15T00:00:00"
const DEFAULT_MISSION_TIME_SEC = "10000.0"
const DEFAULT_N_CLUSTER_TARGETS = "300"

const DEFAULT_BASE_SIGMA_THETA_RAD = "1e-4"
const DEFAULT_BIAS_SWEEP_RAD = "0.0,1e-5,2e-5,5e-5,1e-4"
const DEFAULT_FALSE_ALARM_SWEEP = "0.2,0.4,0.8"
const DEFAULT_MISDETECTION_SWEEP = "0.0,0.02,0.05,0.10"

const RANDOMIZE_SEEDS = false

@inline _truthy(value::AbstractString)::Bool = lowercase(strip(value)) in ("1", "true", "yes", "y", "on")
_string_or_empty(value) = ismissing(value) ? "" : String(value)

function _parse_values(raw::AbstractString)::Vector{Float64}
    return [parse(Float64, strip(value)) for value in split(strip(raw), ",") if !isempty(strip(value))]
end

function _token(value::Float64)::String
    text = replace(@sprintf("%.4g", value), "." => "p")
    text = replace(text, "-" => "m")
    text = replace(text, "+" => "")
    return replace(text, "e" => "e")
end

Base.@kwdef struct SensitivitySpec
    label::Symbol
    parameter::Symbol
    value::Float64
    sigma_theta_rad::Float64
    enable_misdetections::Bool
    misdetection_rate::Float64
    enable_false_alarms::Bool
    false_alarm_rate::Float64
    enable_measurement_bias::Bool
    measurement_bias_rad::Float64
end

function _sensitivity_specs()::Vector{SensitivitySpec}
    base_sigma = parse(Float64, get(ENV, "SPACEAGORA_SWEEP_BASE_SIGMA_THETA_RAD", DEFAULT_BASE_SIGMA_THETA_RAD))
    bias_values = _parse_values(get(ENV, "SPACEAGORA_SWEEP_BIAS_RAD_VALUES", DEFAULT_BIAS_SWEEP_RAD))
    false_alarm_values = _parse_values(get(ENV, "SPACEAGORA_SWEEP_FALSE_ALARM_VALUES", DEFAULT_FALSE_ALARM_SWEEP))
    misdetection_values = _parse_values(get(ENV, "SPACEAGORA_SWEEP_MISDETECTION_VALUES", DEFAULT_MISDETECTION_SWEEP))

    specs = SensitivitySpec[]
    for value in bias_values
        push!(
            specs,
            SensitivitySpec(
                label=Symbol("bias_$(_token(value))"),
                parameter=:bias,
                value=value,
                sigma_theta_rad=base_sigma,
                enable_misdetections=false,
                misdetection_rate=0.0,
                enable_false_alarms=false,
                false_alarm_rate=0.0,
                enable_measurement_bias=value > 0.0,
                measurement_bias_rad=value
            )
        )
    end

    for value in false_alarm_values
        push!(
            specs,
            SensitivitySpec(
                label=Symbol("false_alarm_$(_token(value))"),
                parameter=:false_alarm,
                value=value,
                sigma_theta_rad=base_sigma,
                enable_misdetections=false,
                misdetection_rate=0.0,
                enable_false_alarms=value > 0.0,
                false_alarm_rate=value,
                enable_measurement_bias=false,
                measurement_bias_rad=0.0
            )
        )
    end

    for value in misdetection_values
        push!(
            specs,
            SensitivitySpec(
                label=Symbol("misdetection_$(_token(value))"),
                parameter=:misdetection,
                value=value,
                sigma_theta_rad=base_sigma,
                enable_misdetections=value > 0.0,
                misdetection_rate=value,
                enable_false_alarms=false,
                false_alarm_rate=0.0,
                enable_measurement_bias=false,
                measurement_bias_rad=0.0
            )
        )
    end

    return specs
end

function _sweep_run_config()
    randomize_seeds = _truthy(get(ENV, "SPACEAGORA_RANDOMIZE_SEEDS", string(RANDOMIZE_SEEDS)))
    scenario_seed_int = randomize_seeds ?
                        rand(RandomDevice(), 1:2_000_000_000) :
                        parse(Int, get(ENV, "SPACEAGORA_SCENARIO_SEED", DEFAULT_SCENARIO_SEED))
    sensor_seed_int = randomize_seeds ?
                      scenario_seed_int + 1 :
                      parse(Int, get(ENV, "SPACEAGORA_SENSOR_SEED", string(scenario_seed_int + 1)))
    od_seed_int = randomize_seeds ?
                  scenario_seed_int + 2 :
                  parse(Int, get(ENV, "SPACEAGORA_OBSERVER_OD_SEED", string(scenario_seed_int + 2)))

    return (
        scenario=string(scenario_seed_int),
        sensor=string(sensor_seed_int),
        od=string(od_seed_int),
        start_utc=strip(get(ENV, "SPACEAGORA_START_UTC", DEFAULT_START_UTC)),
        randomize_seeds=randomize_seeds
    )
end

function _case_dir(spec::SensitivitySpec)::String
    return joinpath(OUTPUT_ROOT, String(spec.label))
end

function _metrics_path(spec::SensitivitySpec)::String
    return joinpath(_case_dir(spec), "association_quality_table.csv")
end

function _process_exit_code(process)::Int
    try
        return Int(getproperty(process, :exitcode))
    catch
        return 1
    end
end

function _process_term_signal(process)::Int
    try
        return Int(getproperty(process, :termsignal))
    catch
        return 0
    end
end

function _failure_reason(log_path::String, exit_code::Int, term_signal::Int)::String
    isfile(log_path) || return "run log missing"
    lines = readlines(log_path)
    isempty(lines) && return "run log is empty"

    patterns = (
        "ERROR:",
        "LoadError",
        "BoundsError",
        "DimensionMismatch",
        "OutOfMemory",
        "UndefVarError",
        "MethodError",
        "ArgumentError",
        "AssertionError",
        "Killed",
        "Terminated",
        "signal"
    )
    for line in Iterators.reverse(lines)
        if any(pattern -> occursin(pattern, line), patterns)
            return strip(line)
        end
    end

    signal_hint = term_signal == 9 ? " (SIGKILL: usually external kill or memory pressure)" : ""
    tail_start = max(1, length(lines) - 4)
    tail = strip(join(lines[tail_start:end], " | "))
    return "process exited with code $(exit_code), signal $(term_signal)$(signal_hint). Last log lines: $(tail)"
end

function _run_spec!(spec::SensitivitySpec, run_config; clean_output::Bool)
    case_dir = _case_dir(spec)
    clean_output && rm(case_dir; recursive=true, force=true)
    mkpath(case_dir)
    log_path = joinpath(case_dir, "run.log")
    log_io = open(log_path, "w")
    start_time = time()

    child_env = copy(ENV)
    child_env["SPACEAGORA_COMPARISON_OUTPUT"] = OUTPUT_ROOT
    child_env["SPACEAGORA_NAV_CASE"] = "proposed"
    child_env["SPACEAGORA_NAV_OUTPUT_LABEL"] = String(spec.label)
    child_env["SPACEAGORA_SCENARIO_SEED"] = run_config.scenario
    child_env["SPACEAGORA_SENSOR_SEED"] = run_config.sensor
    child_env["SPACEAGORA_OBSERVER_OD_SEED"] = run_config.od
    child_env["SPACEAGORA_START_UTC"] = run_config.start_utc
    child_env["SPACEAGORA_MISSION_TIME_SEC"] = get(ENV, "SPACEAGORA_MISSION_TIME_SEC", DEFAULT_MISSION_TIME_SEC)
    child_env["SPACEAGORA_N_CLUSTER_TARGETS"] = get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", DEFAULT_N_CLUSTER_TARGETS)
    child_env["SPACEAGORA_SIGMA_THETA_RAD"] = string(spec.sigma_theta_rad)
    child_env["SPACEAGORA_ENABLE_MISDETECTIONS"] = string(spec.enable_misdetections)
    child_env["SPACEAGORA_MISDETECTION_RATE"] = string(spec.misdetection_rate)
    child_env["SPACEAGORA_ENABLE_FALSE_ALARMS"] = string(spec.enable_false_alarms)
    child_env["SPACEAGORA_FALSE_ALARM_RATE"] = string(spec.false_alarm_rate)
    child_env["SPACEAGORA_ENABLE_MEASUREMENT_BIAS"] = string(spec.enable_measurement_bias)
    child_env["SPACEAGORA_MEASUREMENT_BIAS_RAD"] = string(spec.measurement_bias_rad)
    child_env["SPACEAGORA_ENABLE_NAV_TIMING"] = get(ENV, "SPACEAGORA_ENABLE_NAV_TIMING", "false")
    child_env["SPACEAGORA_SAVE_SIMULATION_RESULTS"] = "false"
    child_env["SPACEAGORA_SAVE_TARGET_ESTIMATE_FIELDS"] = "false"
    child_env["SPACEAGORA_SAVE_BUNDLE"] = "0"
    child_env["SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES"] = get(ENV, "SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES", "false")
    child_env["SPACEAGORA_LOG_NAV_EVENTS"] = get(ENV, "SPACEAGORA_LOG_NAV_EVENTS", "false")
    if haskey(ENV, "SPACEAGORA_CASE_JULIA_NUM_THREADS")
        child_env["JULIA_NUM_THREADS"] = ENV["SPACEAGORA_CASE_JULIA_NUM_THREADS"]
    end

    julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
    cmd = setenv(`$julia_exe --project=$REPO_ROOT $SIM_FILE`, child_env)
    println("Running $(spec.label) ($(spec.parameter)=$(spec.value)) -> $(case_dir)")

    process_ok = false
    process_exit_code = 1
    process_term_signal = 0
    failure_reason = ""
    try
        process = run(pipeline(ignorestatus(cmd); stdout=log_io, stderr=log_io))
        process_ok = success(process)
        process_exit_code = _process_exit_code(process)
        process_term_signal = _process_term_signal(process)
    catch err
        println(log_io)
        println(log_io, "Runner error while executing sensitivity case $(spec.label):")
        showerror(log_io, err, catch_backtrace())
        println(log_io)
        process_ok = false
        process_exit_code = 1
    finally
        close(log_io)
    end

    metrics_ok = isfile(_metrics_path(spec))
    ok = process_ok && metrics_ok
    if !process_ok
        failure_reason = _failure_reason(log_path, process_exit_code, process_term_signal)
    elseif !metrics_ok
        failure_reason = "process exited successfully but association_quality_table.csv was not produced"
    end

    elapsed = time() - start_time
    println("Finished $(spec.label): $(ok ? "ok" : "failed") in $(round(elapsed; digits=2)) s")
    ok || println("  reason: $(failure_reason)")

    return (
        String(spec.label),
        String(spec.parameter),
        spec.value,
        spec.sigma_theta_rad,
        spec.enable_misdetections ? 1 : 0,
        spec.misdetection_rate,
        spec.enable_false_alarms ? 1 : 0,
        spec.false_alarm_rate,
        spec.enable_measurement_bias ? 1 : 0,
        spec.measurement_bias_rad,
        run_config.start_utc,
        parse(Int, run_config.scenario),
        parse(Int, run_config.sensor),
        parse(Int, run_config.od),
        ok ? "ok" : "failed",
        process_exit_code,
        process_term_signal,
        elapsed,
        stored_output_path(case_dir, OUTPUT_ROOT),
        stored_output_path(log_path, OUTPUT_ROOT),
        failure_reason
    )
end

function _empty_status_table()::DataFrame
    return DataFrame(
        case=String[],
        parameter=String[],
        value=Float64[],
        sigma_theta_rad=Float64[],
        enable_misdetections=Int[],
        misdetection_rate=Float64[],
        enable_false_alarms=Int[],
        false_alarm_rate=Float64[],
        enable_measurement_bias=Int[],
        measurement_bias_rad=Float64[],
        start_utc=String[],
        scenario_seed=Int[],
        sensor_seed=Int[],
        observer_od_seed=Int[],
        status=String[],
        exit_code=Int[],
        term_signal=Int[],
        elapsed_sec=Float64[],
        results_dir=String[],
        log_path=String[],
        failure_reason=String[]
    )
end

function _run_specs!(specs::Vector{SensitivitySpec}, run_config; clean_output::Bool)::DataFrame
    mkpath(OUTPUT_ROOT)
    live_status_path = joinpath(OUTPUT_ROOT, "sensitivity_run_status_live.csv")
    status_path = joinpath(OUTPUT_ROOT, "sensitivity_run_status.csv")
    skip_existing_ok = _truthy(get(ENV, "SPACEAGORA_SKIP_EXISTING_OK", "false"))
    previous_status = (skip_existing_ok && isfile(status_path)) ? CSV.read(status_path, DataFrame) : DataFrame()
    previous_by_case = Dict{String, DataFrameRow}()
    if nrow(previous_status) > 0 && :case in propertynames(previous_status)
        for row in eachrow(previous_status)
            previous_by_case[String(row.case)] = row
        end
    end

    status_rows = _empty_status_table()
    for spec in specs
        previous_row = get(previous_by_case, String(spec.label), nothing)
        if previous_row !== nothing &&
                String(previous_row.status) == "ok" &&
                isfile(_metrics_path(spec))
            println("Skipping $(spec.label): existing ok result found")
            push!(
                status_rows,
                (
                    String(previous_row.case),
                    String(previous_row.parameter),
                    Float64(previous_row.value),
                    Float64(previous_row.sigma_theta_rad),
                    Int(previous_row.enable_misdetections),
                    Float64(previous_row.misdetection_rate),
                    Int(previous_row.enable_false_alarms),
                    Float64(previous_row.false_alarm_rate),
                    Int(previous_row.enable_measurement_bias),
                    Float64(previous_row.measurement_bias_rad),
                    String(previous_row.start_utc),
                    Int(previous_row.scenario_seed),
                    Int(previous_row.sensor_seed),
                    Int(previous_row.observer_od_seed),
                    String(previous_row.status),
                    Int(previous_row.exit_code),
                    Int(previous_row.term_signal),
                    Float64(previous_row.elapsed_sec),
                    stored_output_path(
                        resolve_output_path(String(previous_row.results_dir), OUTPUT_ROOT),
                        OUTPUT_ROOT
                    ),
                    stored_output_path(
                        resolve_output_path(String(previous_row.log_path), OUTPUT_ROOT),
                        OUTPUT_ROOT
                    ),
                    _string_or_empty(previous_row.failure_reason)
                )
            )
            CSV.write(live_status_path, status_rows)
            continue
        end

        running_row = (
            String(spec.label),
            String(spec.parameter),
            spec.value,
            spec.sigma_theta_rad,
            spec.enable_misdetections ? 1 : 0,
            spec.misdetection_rate,
            spec.enable_false_alarms ? 1 : 0,
            spec.false_alarm_rate,
            spec.enable_measurement_bias ? 1 : 0,
            spec.measurement_bias_rad,
            run_config.start_utc,
            parse(Int, run_config.scenario),
            parse(Int, run_config.sensor),
            parse(Int, run_config.od),
            "running",
            0,
            0,
            0.0,
            stored_output_path(_case_dir(spec), OUTPUT_ROOT),
            stored_output_path(joinpath(_case_dir(spec), "run.log"), OUTPUT_ROOT),
            "case started but no completion status was written"
        )
        push!(status_rows, running_row)
        CSV.write(live_status_path, status_rows)
        status_rows[end, :] = _run_spec!(spec, run_config; clean_output=clean_output)
        CSV.write(live_status_path, status_rows)
    end

    return status_rows
end

function _read_association_table(spec::SensitivitySpec)::DataFrame
    metrics_path = _metrics_path(spec)
    if !isfile(metrics_path)
        return DataFrame(case=String[], parameter=String[], value=Float64[], section=String[], metric=String[], metric_value=Float64[])
    end

    df = CSV.read(metrics_path, DataFrame)
    rename!(df, :value => :metric_value)
    insertcols!(df, 1, :case => fill(String(spec.label), nrow(df)))
    insertcols!(df, 2, :parameter => fill(String(spec.parameter), nrow(df)))
    insertcols!(df, 3, :value => fill(spec.value, nrow(df)))
    return df
end

function _write_sweep_tables(specs::Vector{SensitivitySpec}, status_rows::DataFrame)
    mkpath(OUTPUT_ROOT)
    long_tables = [_read_association_table(spec) for spec in specs]
    long_df = isempty(long_tables) ?
              DataFrame(case=String[], parameter=String[], value=Float64[], section=String[], metric=String[], metric_value=Float64[]) :
              vcat(long_tables...; cols=:union)

    status_path = joinpath(OUTPUT_ROOT, "sensitivity_run_status.csv")
    long_path = joinpath(OUTPUT_ROOT, "sensitivity_metrics_long.csv")

    CSV.write(status_path, status_rows)
    CSV.write(long_path, long_df)

    println()
    println("Sensitivity sweep outputs:")
    println("  status: $(status_path)")
    println("  long metrics: $(long_path)")

    println()
    println("Compact sensitivity summary (one realization per configuration)")
    println("  Cross-run percentile intervals are available only from Monte Carlo analysis.")
    for spec in specs
        case_name = String(spec.label)
        case_metrics = long_df[long_df.case .== case_name, :]
        status_match = status_rows[status_rows.case .== case_name, :status]
        status = isempty(status_match) ? "missing" : String(status_match[1])
        println(
            "  configuration=$(case_name), fault=$(spec.parameter), " *
            "value=$(_compact_number(spec.value)), status=$(status)"
        )

        realized_section, realized_metric, realized_label, realized_unit =
            if spec.parameter === :bias
                (
                    "sensor_errors",
                    "realized_bias_component_rms_rad",
                    "Realized bias RMS",
                    " rad"
                )
            elseif spec.parameter === :misdetection
                (
                    "sensor_errors",
                    "realized_missed_rate_pct",
                    "Realized missed detections",
                    "%"
                )
            else
                (
                    "sensor_errors",
                    "realized_false_alarm_mean_per_observer_epoch",
                    "Realized false alarms / observer-epoch",
                    ""
                )
            end
        realized = _compact_metric(
            case_metrics,
            realized_section,
            realized_metric;
            value_column=:metric_value
        )
        println("    $(realized_label): $(_compact_number(realized))$(realized_unit)")
        _print_compact_navigation_metrics(
            case_metrics;
            value_column=:metric_value
        )

        if spec.parameter === :false_alarm
            fake_tracks = _compact_metric(
                case_metrics,
                "false_alarm",
                "fake_tracks_initialized";
                value_column=:metric_value
            )
            println("    Fake initialized tracks: $(_compact_number(fake_tracks))")
        end
    end
    return (status=status_path, long=long_path)
end

function main()
    specs = _sensitivity_specs()
    run_config = _sweep_run_config()
    clean_output = _truthy(get(ENV, "SPACEAGORA_CLEAN_CASE_OUTPUT", "true"))
    skip_runs = _truthy(get(ENV, "SPACEAGORA_SKIP_RUNS", "false"))

    println("Navigation sensitivity sweep runner")
    println("  time: $(Dates.now(Dates.UTC)) UTC")
    println("  cases: $(length(specs))")
    println("  execution: sequential")
    println("  randomize_seeds: $(run_config.randomize_seeds)")
    println("  start_utc: $(run_config.start_utc)")
    println("  scenario_seed: $(run_config.scenario)")
    println("  sensor_seed: $(run_config.sensor)")
    println("  observer_od_seed: $(run_config.od)")
    println("  mission_time_sec: $(get(ENV, "SPACEAGORA_MISSION_TIME_SEC", DEFAULT_MISSION_TIME_SEC))")
    println("  n_cluster_targets: $(get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", DEFAULT_N_CLUSTER_TARGETS))")
    println("  skip_runs: $(skip_runs)")
    println("  clean_output: $(clean_output)")
    println("  output_root: $(OUTPUT_ROOT)")

    status_rows = if skip_runs
        status_path = joinpath(OUTPUT_ROOT, "sensitivity_run_status.csv")
        isfile(status_path) || error("SPACEAGORA_SKIP_RUNS=true but $(status_path) does not exist")
        CSV.read(status_path, DataFrame)
    else
        _run_specs!(specs, run_config; clean_output=clean_output)
    end
    valid_cases = Set(String(spec.label) for spec in specs)
    status_rows = status_rows[[String(case) in valid_cases for case in status_rows.case], :]

    _write_sweep_tables(specs, status_rows)
    return nothing
end

main()
