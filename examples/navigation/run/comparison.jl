using Dates
using CSV
using DataFrames
using Random

if !isdefined(@__MODULE__, :NavigationPaths)
    include(joinpath(@__DIR__, "..", "paths.jl"))
end
using .NavigationPaths: REPO_ROOT, env_override_path, stored_output_path, transient_navigation_path

include(joinpath(@__DIR__, "compact_summary.jl"))

const SIM_FILE = joinpath(@__DIR__, "..", "..", "navigation.jl")
const OUTPUT_ROOT = env_override_path(
    "SPACEAGORA_COMPARISON_OUTPUT",
    transient_navigation_path("runs", "comparison")
)
# const DEFAULT_SCENARIO_SEED = "1005639413" # ok for 3 obs, 300 tgt, 10000sec
const DEFAULT_SCENARIO_SEED = "1005639417" # ok for alt obs, 300 tgt, 10000sec
const DEFAULT_START_UTC = "2026-06-15T00:00:00"
const DEFAULT_MISSION_TIME_SEC = "10000.0"
const DEFAULT_N_CLUSTER_TARGETS = "300"
const DEFAULT_ENABLE_MISDETECTIONS = "false"
const DEFAULT_ENABLE_FALSE_ALARMS = "false"
const DEFAULT_ENABLE_MEASUREMENT_BIAS = "false"
const DEFAULT_MISDETECTION_RATE = "0.05"
const DEFAULT_FALSE_ALARM_RATE = "0.02"
const DEFAULT_MEASUREMENT_BIAS_RAD = "3e-5"
const DEFAULT_JULIA_HEAP_SIZE_HINT = "4G"
const RANDOMIZE_SEEDS = false
const DEFAULT_BASELINE_THETA_SWEEP = "0.01,0.05"
const DEFAULT_NAV_CASES = (
    :proposed,
    :centralized_oracle,
    :independent_local_da,
    :distributed_oracle_da,
    :baseline_da
)
const VALID_NAV_CASES = (DEFAULT_NAV_CASES..., :no_da)

Base.@kwdef struct CaseSpec
    label::Symbol
    nav_case::Symbol
    baseline_theta_gate_rad::Union{Nothing, Float64} = nothing
end

function _truthy(value::AbstractString)::Bool
    return lowercase(strip(value)) in ("1", "true", "yes", "y", "on")
end

function _baseline_theta_values()::Vector{Float64}
    raw_values = strip(get(ENV, "SPACEAGORA_BASELINE_THETA_SWEEP", DEFAULT_BASELINE_THETA_SWEEP))
    isempty(raw_values) && return Float64[]
    return [parse(Float64, strip(value)) for value in split(raw_values, ",") if !isempty(strip(value))]
end

function _baseline_label(theta_gate_rad::Float64)::Symbol
    token = replace(replace(string(theta_gate_rad), "." => "p"), "-" => "m")
    return Symbol("baseline_da_theta_$(token)")
end

function _case_specs_for(nav_case::Symbol)::Vector{CaseSpec}
    if nav_case === :baseline_da
        theta_values = _baseline_theta_values()
        isempty(theta_values) && return [CaseSpec(label=:baseline_da, nav_case=:baseline_da)]
        return [
            CaseSpec(
                label=_baseline_label(theta),
                nav_case=:baseline_da,
                baseline_theta_gate_rad=theta
            )
            for theta in theta_values
        ]
    end

    return [CaseSpec(label=nav_case, nav_case=nav_case)]
end

function _selected_specs()::Vector{CaseSpec}
    raw_cases = strip(get(ENV, "SPACEAGORA_CASES", ""))
    nav_cases = isempty(raw_cases) ? collect(DEFAULT_NAV_CASES) :
                [Symbol(strip(case)) for case in split(raw_cases, ",") if !isempty(strip(case))]

    baseline_specs_by_label = Dict(spec.label => spec for spec in _case_specs_for(:baseline_da))
    valid = Set(VALID_NAV_CASES)
    invalid = [case for case in nav_cases if !(case in valid) && !haskey(baseline_specs_by_label, case)]
    isempty(invalid) || error(
        "Invalid SPACEAGORA_CASES entries: $(invalid). Valid cases: $(VALID_NAV_CASES) or baseline sweep labels: $(collect(keys(baseline_specs_by_label)))"
    )

    specs = CaseSpec[]
    seen = Set{Symbol}()
    for nav_case in nav_cases
        case_specs = haskey(baseline_specs_by_label, nav_case) ?
                     [baseline_specs_by_label[nav_case]] : _case_specs_for(nav_case)
        for spec in case_specs
            spec.label in seen && continue
            push!(specs, spec)
            push!(seen, spec.label)
        end
    end
    return specs
end

function _case_dir(spec::CaseSpec)::String
    return joinpath(OUTPUT_ROOT, String(spec.label))
end

function _metrics_path(spec::CaseSpec)::String
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

function _signal_hint(term_signal::Int)::String
    term_signal == 9 && return "likely SIGKILL/OOM"
    term_signal == 15 && return "likely SIGTERM/time limit"
    term_signal == 2 && return "interrupted"
    term_signal > 0 && return "terminated by signal $(term_signal)"
    return ""
end

function _exit_code_hint(exit_code::Int, term_signal::Int=0)::String
    signal_hint = _signal_hint(term_signal)
    !isempty(signal_hint) && return signal_hint
    exit_code == 137 && return "likely SIGKILL/OOM"
    exit_code == 143 && return "likely SIGTERM/time limit"
    exit_code == 130 && return "interrupted"
    return ""
end

function _failure_reason(log_path::String, exit_code::Int, term_signal::Int)::String
    isfile(log_path) || return "run log missing"
    lines = readlines(log_path)
    if isempty(lines)
        hint = _exit_code_hint(exit_code, term_signal)
        return isempty(hint) ? "run log is empty" :
               "$(hint); run log is empty (exit_code=$(exit_code), term_signal=$(term_signal))"
    end

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
        "Solve failed",
        "retcode",
        "Killed",
        "Terminated",
        "signal"
    )
    for line in Iterators.reverse(lines)
        if any(pattern -> occursin(pattern, line), patterns)
            return strip(line)
        end
    end

    tail_start = max(1, length(lines) - 4)
    hint = _exit_code_hint(exit_code, term_signal)
    prefix = isempty(hint) ?
             "process exited with code $(exit_code), but log contains no explicit Julia error" :
             "process exited with code $(exit_code) ($(hint)), but log contains no explicit Julia error"
    return "$(prefix). Last log lines: $(strip(join(lines[tail_start:end], " | ")))"
end

function _comparison_run_config()
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
    start_utc = strip(get(ENV, "SPACEAGORA_START_UTC", DEFAULT_START_UTC))

    return (
        scenario=string(scenario_seed_int),
        sensor=string(sensor_seed_int),
        od=string(od_seed_int),
        start_utc=start_utc,
        randomize_seeds=randomize_seeds
    )
end

function _run_case!(spec::CaseSpec, run_config; clean_output::Bool)
    case_dir = _case_dir(spec)
    clean_output && rm(case_dir; recursive=true, force=true)
    mkpath(case_dir)
    log_path = joinpath(case_dir, "run.log")
    log_io = open(log_path, "w")
    start_time = time()

    child_env = copy(ENV)
    child_env["SPACEAGORA_COMPARISON_OUTPUT"] = OUTPUT_ROOT
    child_env["SPACEAGORA_NAV_CASE"] = String(spec.nav_case)
    child_env["SPACEAGORA_NAV_OUTPUT_LABEL"] = String(spec.label)
    if spec.baseline_theta_gate_rad !== nothing
        child_env["SPACEAGORA_BASELINE_ORPHAN_ATTACH_MAX_DTHETA_RAD"] = string(spec.baseline_theta_gate_rad)
    end
    child_env["SPACEAGORA_SCENARIO_SEED"] = run_config.scenario
    child_env["SPACEAGORA_SENSOR_SEED"] = run_config.sensor
    child_env["SPACEAGORA_OBSERVER_OD_SEED"] = run_config.od
    child_env["SPACEAGORA_START_UTC"] = run_config.start_utc
    child_env["SPACEAGORA_MISSION_TIME_SEC"] = get(ENV, "SPACEAGORA_MISSION_TIME_SEC", DEFAULT_MISSION_TIME_SEC)
    child_env["SPACEAGORA_N_CLUSTER_TARGETS"] = get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", DEFAULT_N_CLUSTER_TARGETS)
    child_env["SPACEAGORA_ENABLE_MISDETECTIONS"] = get(ENV, "SPACEAGORA_ENABLE_MISDETECTIONS", DEFAULT_ENABLE_MISDETECTIONS)
    child_env["SPACEAGORA_ENABLE_FALSE_ALARMS"] = get(ENV, "SPACEAGORA_ENABLE_FALSE_ALARMS", DEFAULT_ENABLE_FALSE_ALARMS)
    child_env["SPACEAGORA_ENABLE_MEASUREMENT_BIAS"] = get(ENV, "SPACEAGORA_ENABLE_MEASUREMENT_BIAS", DEFAULT_ENABLE_MEASUREMENT_BIAS)
    child_env["SPACEAGORA_MISDETECTION_RATE"] = get(ENV, "SPACEAGORA_MISDETECTION_RATE", DEFAULT_MISDETECTION_RATE)
    child_env["SPACEAGORA_FALSE_ALARM_RATE"] = get(ENV, "SPACEAGORA_FALSE_ALARM_RATE", DEFAULT_FALSE_ALARM_RATE)
    child_env["SPACEAGORA_MEASUREMENT_BIAS_RAD"] = get(ENV, "SPACEAGORA_MEASUREMENT_BIAS_RAD", DEFAULT_MEASUREMENT_BIAS_RAD)
    child_env["SPACEAGORA_ENABLE_NAV_TIMING"] = spec.nav_case === :proposed ?
        get(ENV, "SPACEAGORA_ENABLE_NAV_TIMING", "true") :
        "false"
    child_env["SPACEAGORA_ENABLE_IOD_ONE_STEP_VALIDATION"] = spec.nav_case in (:proposed, :independent_local_da, :distributed_oracle_da) ?
        get(ENV, "SPACEAGORA_ENABLE_IOD_ONE_STEP_VALIDATION", "true") :
        "false"
    if !_truthy(get(ENV, "SPACEAGORA_COMPARISON_SAVE_TIMESERIES", "false"))
        child_env["SPACEAGORA_SAVE_SIMULATION_RESULTS"] = "false"
        child_env["SPACEAGORA_SAVE_TARGET_ESTIMATE_FIELDS"] = "false"
        child_env["SPACEAGORA_SAVE_BUNDLE"] = "0"
    end
    child_env["SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES"] =
        get(ENV, "SPACEAGORA_SAVE_COMPARISON_DETAILED_TABLES", "false")
    if haskey(ENV, "SPACEAGORA_CASE_JULIA_NUM_THREADS")
        child_env["JULIA_NUM_THREADS"] = ENV["SPACEAGORA_CASE_JULIA_NUM_THREADS"]
    end

    julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
    heap_size_hint = strip(get(ENV, "SPACEAGORA_JULIA_HEAP_SIZE_HINT", DEFAULT_JULIA_HEAP_SIZE_HINT))
    cmd = if isempty(heap_size_hint)
        setenv(`$julia_exe --project=$REPO_ROOT $SIM_FILE`, child_env)
    else
        setenv(`$julia_exe --heap-size-hint=$heap_size_hint --project=$REPO_ROOT $SIM_FILE`, child_env)
    end
    theta_msg = spec.baseline_theta_gate_rad === nothing ? "" :
                " (baseline theta=$(spec.baseline_theta_gate_rad) rad)"
    println("Running $(spec.label) [$(spec.nav_case)]$(theta_msg) -> $(case_dir)")

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
        println(log_io, "Runner error while executing case $(spec.label):")
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
        String(spec.nav_case),
        spec.baseline_theta_gate_rad === nothing ? NaN : spec.baseline_theta_gate_rad,
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

function _run_cases!(specs::Vector{CaseSpec}, run_config; clean_output::Bool)::DataFrame
    mkpath(OUTPUT_ROOT)
    live_status_path = joinpath(OUTPUT_ROOT, "run_status_live.csv")
    status_rows = DataFrame(
        case=String[],
        nav_case=String[],
        baseline_theta_gate_rad=Float64[],
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

    for spec in specs
        running_row = (
            String(spec.label),
            String(spec.nav_case),
            spec.baseline_theta_gate_rad === nothing ? NaN : spec.baseline_theta_gate_rad,
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

        status_rows[end, :] = _run_case!(spec, run_config; clean_output=clean_output)
        CSV.write(live_status_path, status_rows)
    end

    return status_rows
end

function _read_association_table(spec::CaseSpec)::DataFrame
    metrics_path = joinpath(_case_dir(spec), "association_quality_table.csv")
    if !isfile(metrics_path)
        return DataFrame(case=String[], section=String[], metric=String[], value=Float64[])
    end

    df = CSV.read(metrics_path, DataFrame)
    insertcols!(df, 1, :case => fill(String(spec.label), nrow(df)))
    return df
end

function _write_comparison_tables(specs::Vector{CaseSpec}, status_rows::DataFrame)
    mkpath(OUTPUT_ROOT)
    long_tables = [_read_association_table(spec) for spec in specs]
    long_df = isempty(long_tables) ? DataFrame(case=String[], section=String[], metric=String[], value=Float64[]) :
              vcat(long_tables...; cols=:union)

    status_path = joinpath(OUTPUT_ROOT, "run_status.csv")
    long_path = joinpath(OUTPUT_ROOT, "comparison_metrics_long.csv")

    CSV.write(status_path, status_rows)
    CSV.write(long_path, long_df)

    println()
    println("Comparison outputs:")
    println("  status: $(status_path)")
    println("  long metrics: $(long_path)")

    println()
    println("Compact comparison summary (one realization)")
    println("  Cross-run percentile intervals are available only from Monte Carlo analysis.")
    for spec in specs
        case_name = String(spec.label)
        case_metrics = long_df[long_df.case .== case_name, :]
        status_match = status_rows[status_rows.case .== case_name, :status]
        status = isempty(status_match) ? "missing" : String(status_match[1])
        println("  case=$(case_name), status=$(status)")
        _print_compact_navigation_metrics(case_metrics)
    end
    return (status=status_path, long=long_path)
end

function main()
    specs = _selected_specs()
    cases = [spec.label for spec in specs]

    skip_runs = _truthy(get(ENV, "SPACEAGORA_SKIP_RUNS", "false"))
    clean_output = _truthy(get(ENV, "SPACEAGORA_CLEAN_CASE_OUTPUT", "true"))
    run_config = _comparison_run_config()
    mission_time_sec = get(ENV, "SPACEAGORA_MISSION_TIME_SEC", DEFAULT_MISSION_TIME_SEC)
    n_cluster_targets = get(ENV, "SPACEAGORA_N_CLUSTER_TARGETS", DEFAULT_N_CLUSTER_TARGETS)
    enable_misdetections = get(ENV, "SPACEAGORA_ENABLE_MISDETECTIONS", DEFAULT_ENABLE_MISDETECTIONS)
    enable_false_alarms = get(ENV, "SPACEAGORA_ENABLE_FALSE_ALARMS", DEFAULT_ENABLE_FALSE_ALARMS)
    enable_measurement_bias = get(ENV, "SPACEAGORA_ENABLE_MEASUREMENT_BIAS", DEFAULT_ENABLE_MEASUREMENT_BIAS)
    misdetection_rate = get(ENV, "SPACEAGORA_MISDETECTION_RATE", DEFAULT_MISDETECTION_RATE)
    false_alarm_rate = get(ENV, "SPACEAGORA_FALSE_ALARM_RATE", DEFAULT_FALSE_ALARM_RATE)
    measurement_bias_rad = get(ENV, "SPACEAGORA_MEASUREMENT_BIAS_RAD", DEFAULT_MEASUREMENT_BIAS_RAD)
    iod_validation_enabled = get(ENV, "SPACEAGORA_ENABLE_IOD_ONE_STEP_VALIDATION", "true")
    iod_validation_threshold_d2 = get(ENV, "SPACEAGORA_IOD_VALIDATION_MAHAL_MAX_D2", "13.82")
    julia_heap_size_hint = get(ENV, "SPACEAGORA_JULIA_HEAP_SIZE_HINT", DEFAULT_JULIA_HEAP_SIZE_HINT)

    println("Navigation comparison runner")
    println("  time: $(Dates.now(Dates.UTC)) UTC")
    println("  cases: $(cases)")
    println("  baseline_theta_sweep: $(_baseline_theta_values())")
    println("  execution: sequential")
    println("  randomize_seeds: $(run_config.randomize_seeds)")
    println("  start_utc: $(run_config.start_utc)")
    println("  scenario_seed: $(run_config.scenario)")
    println("  sensor_seed: $(run_config.sensor)")
    println("  observer_od_seed: $(run_config.od)")
    println("  mission_time_sec: $(mission_time_sec)")
    println("  n_cluster_targets: $(n_cluster_targets)")
    println("  sensor_errors:")
    println("    missed_detections=$(enable_misdetections), rate=$(misdetection_rate)")
    println("    false_alarms=$(enable_false_alarms), rate=$(false_alarm_rate)")
    println("    measurement_bias=$(enable_measurement_bias), bias_rad=$(measurement_bias_rad)")
    println("  iod_one_step_validation: proposed_only=true, proposed_enabled=$(_truthy(iod_validation_enabled)), threshold_d2=$(iod_validation_threshold_d2)")
    println("  julia_heap_size_hint: $(isempty(strip(julia_heap_size_hint)) ? "disabled" : julia_heap_size_hint)")
    println("  skip_runs: $(skip_runs)")
    println("  clean_output: $(clean_output)")
    println("  output_root: $(OUTPUT_ROOT)")

    status_rows = if skip_runs
        DataFrame(
            case=String.(cases),
            nav_case=String.([spec.nav_case for spec in specs]),
            baseline_theta_gate_rad=[
                spec.baseline_theta_gate_rad === nothing ? NaN : spec.baseline_theta_gate_rad
                for spec in specs
            ],
            start_utc=fill(run_config.start_utc, length(cases)),
            scenario_seed=fill(parse(Int, run_config.scenario), length(cases)),
            sensor_seed=fill(parse(Int, run_config.sensor), length(cases)),
            observer_od_seed=fill(parse(Int, run_config.od), length(cases)),
            status=fill("skipped", length(cases)),
            exit_code=fill(0, length(cases)),
            term_signal=fill(0, length(cases)),
            elapsed_sec=fill(0.0, length(cases)),
            results_dir=[
                stored_output_path(_case_dir(spec), OUTPUT_ROOT) for spec in specs
            ],
            log_path=[
                stored_output_path(joinpath(_case_dir(spec), "run.log"), OUTPUT_ROOT)
                for spec in specs
            ],
            failure_reason=fill("", length(cases))
        )
    else
        _run_cases!(specs, run_config; clean_output=clean_output)
    end

    _write_comparison_tables(specs, status_rows)
    return nothing
end

main()
