const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance")
const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH_HARMONICS_FILE = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
const EARTH_GRAM_SURROGATE_FILE = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "Earth", "earth_surrogate.jls")
const PERF_BASELINE_SCENARIO = "single_j2"
const _PERF_POLICY_ENV_NAMES = (
    "SPACEAGORA_OUTER_PARALLEL_ACTIVE",
    "SPACEAGORA_DENSITY_CALLBACK_PARALLEL",
    "SPACEAGORA_CONTROL_CALLBACK_PARALLEL",
    "SPACEAGORA_MULTIBODY_PARALLEL",
)
const _PERF_POLICY_ENV_BASELINE = Dict{String, Union{Nothing, String}}(
    name => (haskey(ENV, name) ? String(ENV[name]) : nothing)
    for name in _PERF_POLICY_ENV_NAMES
)
const _PERF_THREADS_BACKEND_WARNING_EMITTED = Ref(false)

include(joinpath(REPO_ROOT, "src", "parallel", "ParallelProfiles.jl"))
using .ParallelProfiles

using CSV
using DataFrames
using Dates
using Distributed
using LinearAlgebra
using Random
using Sockets
using SPICE
using StaticArrays
using Statistics

if myid() == 1
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
    using Plots
end

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))
using .SimulationModel

# run_simulation.jl expects quat_mult in the including scope.
const quat_mult = SimulationModel.quat_mult
include(joinpath(REPO_ROOT, "src", "simulation", "run_simulation.jl"))

Base.@kwdef struct ProfileSpec
    name::String
    repeats::Int
    warmup::Int
    max_attempts::Int
    mission_short_s::Float64
    mission_long_s::Float64
    montecarlo_samples::Int
    montecarlo_mission_s::Float64
end

Base.@kwdef struct BenchmarkCase
    name::String
    category::String
    description::String
    args_template::SimulationConfiguration
    run_in_quick::Bool = true
    solver_mode_override::Union{Nothing, String} = nothing
    split_imex_solver_override::Union{Nothing, String} = nothing
    entry_target_count_override::Union{Nothing, Int} = nothing
    env_overrides::Vector{Pair{String, String}} = Pair{String, String}[]
end

@inline _safe_div(num::Float64, den::Float64) = den > 0.0 ? num / den : NaN

@inline function _alloc_calls(gcstats::Base.GC_Diff)::Int
    return gcstats.malloc + gcstats.realloc + gcstats.poolalloc + gcstats.bigalloc
end

@inline function _destat_int(sol, field::Symbol)
    if hasproperty(sol.destats, field)
        value = getproperty(sol.destats, field)
        return value isa Integer ? Int(value) : try
            Int(value)
        catch
            missing
        end
    end
    return missing
end

@inline function _solve_success(sol)::Bool
    retcode = string(sol.retcode)
    # For paper/runtime benchmarking, only full mission completion is valid.
    # "Terminated" usually means an early-stop callback path (e.g., impact),
    # which makes timing samples non-comparable.
    return retcode == "Success"
end

@inline function _solve_success_for_case(sol, case::BenchmarkCase)::Bool
    if _solve_success(sol)
        return true
    end
    if string(sol.retcode) == "Terminated" &&
       case.category == "entry" &&
       !(isnothing(case.entry_target_count_override)) &&
       case.entry_target_count_override > 0
        # Entry cases intentionally terminate after target interface crossings.
        return true
    end
    return false
end

@inline function _effector_signature(effectors::Tuple)
    isempty(effectors) && return "none"
    return join([string(nameof(typeof(e))) for e in effectors], "+")
end

@inline function _safe_unique_join(values; delimiter::String=",")
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return join(sort(unique(string.(vec))), delimiter)
end

function _primary_terminal_state_metrics(sol)
    pos_norm_m = missing
    vel_norm_mps = missing
    mass_kg = missing

    if length(sol.u) == 0
        return (pos_norm_m=pos_norm_m, vel_norm_mps=vel_norm_mps, mass_kg=mass_kg)
    end
    final_state = sol.u[end]
    if !(hasproperty(final_state, :sc))
        return (pos_norm_m=pos_norm_m, vel_norm_mps=vel_norm_mps, mass_kg=mass_kg)
    end
    states = final_state.sc
    isempty(states) && return (pos_norm_m=pos_norm_m, vel_norm_mps=vel_norm_mps, mass_kg=mass_kg)
    primary = states[1]

    if hasproperty(primary, :pos)
        value = norm(primary.pos)
        if isfinite(value)
            pos_norm_m = Float64(value)
        end
    end
    if hasproperty(primary, :vel)
        value = norm(primary.vel)
        if isfinite(value)
            vel_norm_mps = Float64(value)
        end
    end
    if hasproperty(primary, :mass)
        value = Float64(primary.mass)
        if isfinite(value)
            mass_kg = value
        end
    end

    return (pos_norm_m=pos_norm_m, vel_norm_mps=vel_norm_mps, mass_kg=mass_kg)
end

@inline _perf_error_text(err) = sprint(showerror, err)

@inline function _perf_strict_errors()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_STRICT_ERRORS", "0")))
    return raw in ("1", "true", "yes", "on")
end

@inline function _perf_warmup_logs()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_WARMUP_LOGS", "0")))
    return raw in ("1", "true", "yes", "on")
end

@inline function _perf_solver_mode_env()::String
    return _perf_solver_mode_env("quick")
end

@inline function _perf_default_solver_mode(profile_name::String)::String
    normalized = lowercase(strip(profile_name))
    # Keep benchmark defaults close to production behavior for both profiles.
    # Individual cases can still pin solver_mode_override when needed.
    if normalized == "full"
        return "auto_stiff"
    end
    # Quick profile is used for threshold tuning and should stay near production defaults.
    return "auto_stiff"
end

@inline function _perf_solver_mode_env(profile_name::String)::String
    mode = strip(get(ENV, "SPACEAGORA_PERF_SOLVER_MODE", ""))
    if isempty(mode)
        mode = strip(get(ENV, "SPACEAGORA_SOLVER_MODE", ""))
    end
    if isempty(mode)
        mode = _perf_default_solver_mode(profile_name)
    end
    return mode
end

@inline function _perf_stream_orbit_logs()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_STREAM_ORBIT_LOGS", "1")))
    return raw in ("1", "true", "yes", "on")
end

function _run_perf_simulation(
    args_run;
    return_solution::Bool,
    return_solver_metadata::Bool=false,
    profile_name::String="quick",
    solver_mode_override::Union{Nothing, String}=nothing,
    split_imex_solver_override::Union{Nothing, String}=nothing,
    entry_target_count_override::Union{Nothing, Int}=nothing
)
    solver_mode = isnothing(solver_mode_override) ? _perf_solver_mode_env(profile_name) : solver_mode_override
    split_solver_env = isnothing(split_imex_solver_override) ? nothing : String(split_imex_solver_override)
    entry_target_env = if isnothing(entry_target_count_override)
        nothing
    else
        value = Int(entry_target_count_override)
        value >= 0 || throw(ArgumentError("entry_target_count_override must be >= 0, got $value"))
        string(value)
    end
    return withenv(
        "SPACEAGORA_SOLVER_MODE" => solver_mode,
        "SPACEAGORA_SPLIT_IMEX_SOLVER" => split_solver_env,
        "SPACEAGORA_ENTRY_TARGET_COUNT" => entry_target_env
    ) do
        run_simulation(
            args_run;
            isolate_state=false,
            return_solution=return_solution,
            return_solver_metadata=return_solver_metadata
        )
    end
end

@inline function _perf_stack_head(bt)::String
    bt === nothing && return "stack=unavailable"
    st = stacktrace(bt)
    if isempty(st)
        return "stack=empty"
    end
    for frame in st
        file = String(frame.file)
        if !(occursin("/julia/", file) || occursin("boot.jl", file) || occursin("task.jl", file))
            return string(file, ":", frame.line, " in ", frame.func)
        end
    end
    frame = st[1]
    return string(String(frame.file), ":", frame.line, " in ", frame.func)
end

@inline function _solve_retcode_from_error(err)::Union{Missing, String}
    msg = _perf_error_text(err)
    m = match(r"retcode=([A-Za-z0-9_]+)", msg)
    return m === nothing ? missing : String(m.captures[1])
end

@inline function _profile_from_name(name::String)::ProfileSpec
    if name == "full"
        return ProfileSpec(
            name="full",
            repeats=5,
            warmup=2,
            max_attempts=4,
            mission_short_s=3600.0,
            mission_long_s=14400.0,
            montecarlo_samples=50,
            montecarlo_mission_s=3600.0
        )
    elseif name == "quick"
        return ProfileSpec(
            name="quick",
            repeats=3,
            warmup=1,
            max_attempts=4,
            mission_short_s=1800.0,
            mission_long_s=7200.0,
            montecarlo_samples=50,
            montecarlo_mission_s=1800.0
        )
    else
        throw(ArgumentError("Unsupported profile '$name'. Valid values: quick, full."))
    end
end

function parse_cli()
    profile_name = lowercase(get(ENV, "SPACEAGORA_PERF_PROFILE", "quick"))
    outdir = get(ENV, "SPACEAGORA_PERF_OUTDIR", DEFAULT_OUTPUT_DIR)

    for arg in ARGS
        if arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        else
            throw(ArgumentError("Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..."))
        end
    end

    return _profile_from_name(profile_name), abspath(outdir)
end

@inline function perf_parallel_enabled()::Bool
    mode = lowercase(strip(get(ENV, "SPACEAGORA_PERF_PARALLEL", "auto")))
    return !(mode in ("0", "false", "off", "no"))
end

@inline perf_threads_available()::Bool = Threads.nthreads() > 1


@inline function perf_parallel_backend()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_PERF_PARALLEL_BACKEND", "auto")))
    if mode in ("none", "serial", "off", "0", "false", "no")
        return :none
    elseif mode in ("threads", "thread")
        if !perf_parallel_enabled()
            return :none
        end
        if perf_threads_available()
            return :threads
        end
        if !_PERF_THREADS_BACKEND_WARNING_EMITTED[]
            @warn "[perf] SPACEAGORA_PERF_PARALLEL_BACKEND=threads requested but JULIA_NUM_THREADS=$(Threads.nthreads()). Falling back to serial backend (:none). Set JULIA_NUM_THREADS>1 or use SPACEAGORA_PERF_PARALLEL_BACKEND=auto/process."
            _PERF_THREADS_BACKEND_WARNING_EMITTED[] = true
        end
        return :none
    elseif mode in ("process", "processes", "distributed", "pmap")
        return :process
    elseif mode == "auto"
        return :auto
    else
        throw(ArgumentError("Unsupported SPACEAGORA_PERF_PARALLEL_BACKEND='$mode'. Use one of: threads, process, none, auto."))
    end
end

@inline _threads_or_none_backend()::Symbol = perf_threads_available() ? :threads : :none

Base.@kwdef struct ParallelPriorityPlan
    outer_route::Symbol = :none
    density_mode::String = "off"
    control_mode::String = "off"
    multibody_mode::String = "off"
end

const _PERF_OUTER_ROUTE_STATE = Ref(ParallelProfiles.OuterRouteState())

@inline function _parse_positive_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    return max(1, value)
end

@inline function _parse_nonnegative_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    return max(0, value)
end

@inline function _parse_nonnegative_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a number, got '$raw'"))
    end
    return max(0.0, value)
end

@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _priority_inner_sat_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_INNER_SAT_THRESHOLD", 8)
end

@inline function _priority_inner_link_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_INNER_LINK_THRESHOLD", 12)
end

@inline function _priority_outer_light_sat_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_OUTER_LIGHT_SAT_THRESHOLD", 2)
end

@inline function _priority_outer_light_link_threshold()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_OUTER_LIGHT_LINK_THRESHOLD", 4)
end

@inline function _priority_outer_light_mission_threshold_s()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_PRIORITY_OUTER_LIGHT_MISSION_THRESHOLD_S", 14_400.0)
end

@inline function _priority_spice_constellation_process_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_PRIORITY_SPICE_CONSTELLATION_PROCESS", true)
end

@inline function _priority_spice_constellation_min_sats()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_PRIORITY_SPICE_CONSTELLATION_MIN_SATS", 4)
end

@inline function _outer_route_adaptive_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE", true)
end

@inline function _outer_route_min_samples()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_OUTER_ROUTE_MIN_SAMPLES", 2)
end

@inline function _outer_route_failure_penalty_s()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_OUTER_ROUTE_FAILURE_PENALTY_S", 120.0)
end

@inline function _outer_route_mc_process_min_samples()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_MC_PROCESS_MIN_SAMPLES", 16)
end

@inline function _outer_route_mc_process_min_mission_s()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MC_PROCESS_MIN_MISSION_S", 3600.0)
end

@inline function _outer_route_trace_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_OUTER_ROUTE_TRACE", false)
end

@inline function _outer_route_state_persist_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_OUTER_ROUTE_STATE_PERSIST", true)
end

@inline function _outer_route_state_reset_requested()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_OUTER_ROUTE_STATE_RESET", false)
end

@inline function _control_stress_repeats_full()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL", 3)
end

@inline function _control_stress_warmup_full()::Int
    return _parse_nonnegative_int_env("SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL", 1)
end

@inline function _include_control_stress_per_orbit()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_INCLUDE_CONTROL_STRESS_PER_ORBIT", true)
end

@inline function _split_rollout_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_SPLIT_ROLLOUT_GATE", false)
end

@inline function _split_rollout_enforce()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_SPLIT_ROLLOUT_ENFORCE", false)
end

@inline function _split_rollout_case_names()::Vector{String}
    raw = strip(get(
        ENV,
        "SPACEAGORA_PERF_SPLIT_ROLLOUT_CASES",
        "single_orientation_aero,proximity_2sat_orientation_fullstack_gnc_highrate,multi_4_gravity"
    ))
    tokens = String[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(tokens, t)
    end
    return unique(tokens)
end

@inline function _split_rollout_solver_variants()::Vector{String}
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_SPLIT_ROLLOUT_SOLVERS", "kencarp4,kencarp47")))
    tokens = String[]
    for token in split(raw, ",")
        t = lowercase(strip(token))
        isempty(t) && continue
        if !(t in ("kencarp4", "kencarp47", "kencarp58"))
            throw(ArgumentError(
                "SPACEAGORA_PERF_SPLIT_ROLLOUT_SOLVERS token '$t' is unsupported. " *
                "Use comma-separated values from: kencarp4, kencarp47, kencarp58."
            ))
        end
        push!(tokens, t)
    end
    isempty(tokens) && return ["kencarp4"]
    return unique(tokens)
end

@inline function _split_rollout_max_slowdown_ratio()::Float64
    return max(eps(Float64), _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_MAX_SLOWDOWN_RATIO", 1.15))
end

@inline function _split_rollout_pos_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_POS_REL_TOL", 5e-4)
end

@inline function _split_rollout_vel_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_VEL_REL_TOL", 5e-4)
end

@inline function _split_rollout_q_angle_tol_rad()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_Q_ANGLE_TOL_RAD", 5e-4)
end

@inline function _split_rollout_omega_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_SPLIT_OMEGA_REL_TOL", 1e-3)
end

@inline function _split_rollout_sample_count()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_SPLIT_TRAJ_SAMPLES", 128)
end

@inline function _multirate_rollout_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_MULTIRATE_ROLLOUT_GATE", false)
end

@inline function _multirate_rollout_enforce()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_MULTIRATE_ROLLOUT_ENFORCE", false)
end

@inline function _multirate_rollout_case_names()::Vector{String}
    raw = strip(get(
        ENV,
        "SPACEAGORA_PERF_MULTIRATE_ROLLOUT_CASES",
        "single_orientation_aero,proximity_2sat_orientation_fullstack_gnc_highrate,multi_4_gravity"
    ))
    tokens = String[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(tokens, t)
    end
    return unique(tokens)
end

@inline function _multirate_rollout_max_slowdown_ratio()::Float64
    return max(eps(Float64), _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_MAX_SLOWDOWN_RATIO", 1.25))
end

@inline function _multirate_rollout_pos_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_POS_REL_TOL", 7.5e-4)
end

@inline function _multirate_rollout_vel_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_VEL_REL_TOL", 7.5e-4)
end

@inline function _multirate_rollout_q_angle_tol_rad()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_Q_ANGLE_TOL_RAD", 7.5e-4)
end

@inline function _multirate_rollout_omega_rel_tol()::Float64
    return _parse_nonnegative_float_env("SPACEAGORA_PERF_MULTIRATE_OMEGA_REL_TOL", 1.5e-3)
end

@inline function _multirate_rollout_sample_count()::Int
    return _parse_positive_int_env("SPACEAGORA_PERF_MULTIRATE_TRAJ_SAMPLES", 128)
end

@inline function _multirate_env_float_or_missing(name::String)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return missing
    parsed = try
        parse(Float64, raw)
    catch
        return missing
    end
    return parsed
end

@inline function _multirate_env_int_or_missing(name::String)
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return missing
    parsed = try
        parse(Int, raw)
    catch
        return missing
    end
    return parsed
end

@inline function _multirate_rollout_setting_snapshot()
    slow_solver = lowercase(strip(get(ENV, "SPACEAGORA_MULTIRATE_SLOW_SOLVER", "tsit5")))
    fast_solver = lowercase(strip(get(ENV, "SPACEAGORA_MULTIRATE_FAST_SOLVER", "auto_stiff")))
    slow_dt_s = _multirate_env_float_or_missing("SPACEAGORA_MULTIRATE_SLOW_DT_S")
    fast_substeps = _multirate_env_int_or_missing("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS")
    return (
        slow_solver=slow_solver,
        fast_solver=fast_solver,
        slow_dt_s=slow_dt_s,
        fast_substeps=fast_substeps
    )
end

@inline function _machine_parallel_class()::Symbol
    cpu_threads = Sys.CPU_THREADS
    if cpu_threads >= 24
        return :large
    elseif cpu_threads >= 12
        return :medium
    end
    return :small
end

@inline function _hardware_class_name()::String
    override_raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_HARDWARE_CLASS", "auto")))
    if override_raw == "auto" || isempty(override_raw)
        return String(_machine_parallel_class())
    end
    return override_raw
end

@inline function _machine_label()::String
    raw = strip(get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", ""))
    return isempty(raw) ? gethostname() : raw
end

@inline function _active_parallel_profile_token()::String
    raw = strip(get(ENV, "SPACEAGORA_PARALLEL_PROFILE", "default"))
    return isempty(raw) ? "default" : raw
end

@inline function _active_parallel_machine_token()::String
    raw = strip(get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", "default"))
    return isempty(raw) ? "default" : raw
end

@inline function _safe_path_token(raw::AbstractString)::String
    token = lowercase(strip(String(raw)))
    token = replace(token, r"[^a-z0-9._-]+" => "_")
    return isempty(token) ? "default" : token
end

@inline function _outer_route_state_cache_path(spec::ProfileSpec, outdir::String)::String
    override_raw = strip(get(ENV, "SPACEAGORA_PERF_OUTER_ROUTE_STATE_PATH", ""))
    if !isempty(override_raw)
        return normpath(isabspath(override_raw) ? override_raw : joinpath(REPO_ROOT, override_raw))
    end
    profile_token = _safe_path_token(spec.name)
    machine_token = _safe_path_token(_machine_label())
    hardware_token = _safe_path_token(_hardware_class_name())
    return joinpath(
        outdir,
        "outer_route_state",
        "outer_route_state_$(profile_token)_$(hardware_token)_$(machine_token).toml"
    )
end

@inline function _cpu_name_string()::String
    return hasproperty(Sys, :CPU_NAME) ? String(getproperty(Sys, :CPU_NAME)) : "unknown"
end

@inline function _runtime_hardware_snapshot()
    return (
        hardware_class=_hardware_class_name(),
        machine_label=_machine_label(),
        host_name=gethostname(),
        cpu_name=_cpu_name_string(),
        cpu_threads=Sys.CPU_THREADS,
        julia_threads=Threads.nthreads(),
        os=string(Sys.KERNEL),
        arch=string(Sys.ARCH)
    )
end

function _inner_hint_layer_report_df(spec::ProfileSpec, hw)::DataFrame
    rows = SimulationModel.ParallelPolicy.hint_layer_stats_snapshot(
        profile=_active_parallel_profile_token(),
        machine=_active_parallel_machine_token()
    )
    stamp_utc = string(now(UTC))
    if isempty(rows)
        return DataFrame(
            runtime_profile=String[],
            parallel_profile=String[],
            hint_machine=String[],
            hardware_class=String[],
            machine_label=String[],
            host_name=String[],
            layer=String[],
            signature_count=Int[],
            choice_count=Int[],
            samples_total=Int[],
            successes_total=Int[],
            failures_total=Int[],
            elapsed_mean_ns=Float64[],
            elapsed_std_ns=Float64[],
            confidence_mean=Float64[],
            regret_mean_ns=Float64[],
            exploration_c=Float64[],
            min_samples=Int[],
            state_path=String[],
            timestamp_utc=String[]
        )
    end
    return DataFrame([
        (
            runtime_profile=spec.name,
            parallel_profile=String(row.profile),
            hint_machine=String(row.machine),
            hardware_class=hw.hardware_class,
            machine_label=hw.machine_label,
            host_name=hw.host_name,
            layer=String(row.layer),
            signature_count=Int(row.signature_count),
            choice_count=Int(row.choice_count),
            samples_total=Int(row.samples_total),
            successes_total=Int(row.successes_total),
            failures_total=Int(row.failures_total),
            elapsed_mean_ns=Float64(row.elapsed_mean_ns),
            elapsed_std_ns=Float64(row.elapsed_std_ns),
            confidence_mean=Float64(row.confidence_mean),
            regret_mean_ns=Float64(row.regret_mean_ns),
            exploration_c=Float64(row.exploration_c),
            min_samples=Int(row.min_samples),
            state_path=String(row.state_path),
            timestamp_utc=stamp_utc
        )
        for row in rows
    ])
end

@inline function _total_links(case::BenchmarkCase)::Int
    total = 0
    @inbounds for sc in case.args_template.dynamics_model.spacecraft
        total += length(sc.links)
    end
    return total
end

@inline function _has_atmo_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
            return true
        end
    end
    return false
end

@inline function _has_nbody_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa NBodyGravityModel
            return true
        end
    end
    return false
end

@inline function _has_srp_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa SolarRadiationPressureModel && effector.A > 0.0
            return true
        end
    end
    return false
end

@inline function _max_harmonics_degree(effectors::Tuple)::Int
    degree = 0
    @inbounds for effector in effectors
        if effector isa GravitationalHarmonicsModel
            degree = max(degree, effector.L)
        end
    end
    return degree
end

@inline function _max_links_per_sat(case::BenchmarkCase)::Int
    max_links = 1
    @inbounds for sc in case.args_template.dynamics_model.spacecraft
        max_links = max(max_links, length(sc.links))
    end
    return max_links
end

@inline function _density_model_family(density_model)::String
    if density_model isa NoAtmosphereModel
        return "none"
    elseif density_model isa GRAMAtmosphereModelSurrogate
        return "gram_surrogate"
    elseif density_model isa GRAMAtmosphereModel
        return "gram_point"
    elseif density_model isa ExponentialAtmosphereModel
        return "exponential"
    elseif density_model isa PolynomialFitAtmosphereModel
        return "polyfit"
    elseif density_model isa NRLMSISE00AtmosphereModel
        return "nrlmsise00"
    end
    return lowercase(string(nameof(typeof(density_model))))
end

@inline function _read_bool_property_safe(obj, name::Symbol)::Union{Nothing, Bool}
    if !hasproperty(obj, name)
        return nothing
    end
    value = try
        getproperty(obj, name)
    catch
        return nothing
    end
    return value isa Bool ? value : nothing
end

@inline function _gram_surrogate_flag(density_model)::Bool
    if density_model isa GRAMAtmosphereModelSurrogate
        return true
    end
    # Fallback for custom wrappers exposing a boolean surrogate toggle.
    flag = _read_bool_property_safe(density_model, :gram_offline_surrogate)
    return isnothing(flag) ? false : flag
end

@inline function _gram_static_grid_flag(density_model)::Bool
    if !(density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate)
        return false
    end
    for key in (:use_static_grid, :static_grid, :gram_static_grid)
        flag = _read_bool_property_safe(density_model, key)
        if !(flag === nothing)
            return flag
        end
    end
    return false
end

@inline function _env_bool_token(raw::AbstractString)::Union{Nothing, Bool}
    token = lowercase(strip(String(raw)))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    return nothing
end

@inline function _case_env_value(case::BenchmarkCase, name::String)::Union{Nothing, String}
    for pair in case.env_overrides
        if first(pair) == name
            return String(last(pair))
        end
    end
    return nothing
end

@inline function _env_bool_override_for_case(case::BenchmarkCase, name::String)::Union{Nothing, Bool}
    value = _case_env_value(case, name)
    value === nothing && return nothing
    return _env_bool_token(value)
end

@inline function _gram_static_grid_flag_for_case(case::BenchmarkCase, density_model)::Bool
    if !(density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate)
        return false
    end
    override = _env_bool_override_for_case(case, "SPACEAGORA_GRAM_STATIC_GRID")
    if !(override === nothing)
        return override
    end
    return _gram_static_grid_flag(density_model)
end

@inline function _gram_track_cache_mode_for_case(case::BenchmarkCase)::String
    override = _case_env_value(case, "SPACEAGORA_GRAM_TRACK_CACHE")
    if !(override === nothing)
        token = lowercase(strip(override))
        if token in ("on", "off", "auto")
            return token
        end
    end
    return lowercase(strip(get(ENV, "SPACEAGORA_GRAM_TRACK_CACHE", "off")))
end

@inline function _density_backend_bucket(
    density_family::AbstractString,
    gram_surrogate_enabled::Bool,
    gram_static_grid_enabled::Bool,
    gram_track_cache_mode::AbstractString
)::String
    family = lowercase(strip(String(density_family)))
    cache_on = lowercase(strip(String(gram_track_cache_mode))) in ("on", "auto")
    if family in ("gram_point", "gram")
        return gram_static_grid_enabled ? "gram_static_grid_or_cached_surrogate" : "gram_point_to_point"
    elseif family in ("gram_surrogate", "gram_offline_surrogate")
        if gram_static_grid_enabled || cache_on
            return "gram_static_grid_or_cached_surrogate"
        end
        return gram_surrogate_enabled ? "gram_surrogate" : "gram_point_to_point"
    elseif family == "none"
        return "non_gram"
    elseif occursin("gram", family)
        return gram_static_grid_enabled ? "gram_static_grid_or_cached_surrogate" : "gram_point_to_point"
    end
    return "non_gram"
end

@inline function _solver_mode_for_outer_route(case::BenchmarkCase, profile_name::String)::String
    if !(case.solver_mode_override === nothing)
        return String(case.solver_mode_override)
    end
    if !(case.split_imex_solver_override === nothing)
        return "split_imex:" * String(case.split_imex_solver_override)
    end
    return _perf_solver_mode_env(profile_name)
end

@inline function _min_positive_rate(rates::AbstractVector{<:Real})::Float64
    best = Inf
    @inbounds for rate in rates
        value = Float64(rate)
        if isfinite(value) && value > 0.0
            best = min(best, value)
        end
    end
    return isfinite(best) ? best : 0.0
end

@inline function _has_aero_dynamic_effector(effectors::Tuple)::Bool
    return _has_atmo_dynamic_effector(effectors)
end

@inline function _dynamic_effector_cost_class(effectors::Tuple, control_effector_count::Int)::String
    effector_count = length(effectors)
    has_nbody = _has_nbody_dynamic_effector(effectors)
    has_srp = _has_srp_dynamic_effector(effectors)
    has_aero = _has_aero_dynamic_effector(effectors)
    harmonics_degree = _max_harmonics_degree(effectors)

    if has_nbody || harmonics_degree >= 20 || effector_count >= 5
        return "heavy"
    end
    if has_srp || has_aero || harmonics_degree > 0 || effector_count >= 2 || control_effector_count >= 2
        return "medium"
    end
    return "light"
end

@inline function _case_outer_threads_safe(case::BenchmarkCase)::Bool
    # Guard against native-library aborts seen under outer-threaded case concurrency.
    density_model = case.args_template.environment_model.density_model
    if density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate
        return false
    end
    return !_has_nbody_dynamic_effector(case.args_template.dynamics_model.dynamic_effectors)
end

function _split_threadsafe_indices(
    selected::Vector{BenchmarkCase},
    indices::Vector{Int}
)::Tuple{Vector{Int}, Vector{Int}}
    thread_indices = Int[]
    serial_indices = Int[]
    for idx in indices
        if _case_outer_threads_safe(selected[idx])
            push!(thread_indices, idx)
        else
            push!(serial_indices, idx)
        end
    end
    return thread_indices, serial_indices
end

function _reset_outer_route_history!()
    _PERF_OUTER_ROUTE_STATE[] = ParallelProfiles.OuterRouteState()
    return nothing
end

function _load_outer_route_history!(spec::ProfileSpec, outdir::String)
    _reset_outer_route_history!()
    persist = _outer_route_state_persist_enabled() && _outer_route_adaptive_enabled()
    path = _outer_route_state_cache_path(spec, outdir)
    if !persist
        return (persist=false, reset_requested=false, path=path, loaded_rows=0, loaded_signatures=0)
    end
    reset_requested = _outer_route_state_reset_requested()
    if reset_requested
        return (persist=true, reset_requested=true, path=path, loaded_rows=0, loaded_signatures=0)
    end
    loaded = try
        ParallelProfiles.load_outer_route_state!(_PERF_OUTER_ROUTE_STATE[], path; replace=true)
    catch err
        @warn "[perf] Failed to load outer-route adaptive cache; starting with empty history." path=path exception=(err, catch_backtrace())
        (path=path, rows=0, signatures=0)
    end
    return (
        persist=true,
        reset_requested=false,
        path=loaded.path,
        loaded_rows=loaded.rows,
        loaded_signatures=loaded.signatures
    )
end

function _save_outer_route_history!(spec::ProfileSpec, path::String, enabled::Bool)
    enabled || return (path=path, rows=0, signatures=0)
    metadata = Dict{String, Any}(
        "profile" => spec.name,
        "machine_label" => _machine_label(),
        "hardware_class" => _hardware_class_name(),
        "cpu_threads" => Sys.CPU_THREADS,
        "julia_threads" => Threads.nthreads()
    )
    saved = try
        ParallelProfiles.save_outer_route_state(
            _PERF_OUTER_ROUTE_STATE[],
            path;
            metadata=metadata
        )
    catch err
        @warn "[perf] Failed to persist outer-route adaptive cache." path=path exception=(err, catch_backtrace())
        return (path=path, rows=0, signatures=0)
    end
    return saved
end

@inline function _outer_route_tuning()::ParallelProfiles.OuterRouteTuning
    return ParallelProfiles.OuterRouteTuning(
        inner_sat_threshold=_priority_inner_sat_threshold(),
        inner_link_threshold=_priority_inner_link_threshold(),
        outer_light_sat_threshold=_priority_outer_light_sat_threshold(),
        outer_light_link_threshold=_priority_outer_light_link_threshold(),
        outer_light_mission_threshold_s=_priority_outer_light_mission_threshold_s(),
        spice_constellation_process_enabled=_priority_spice_constellation_process_enabled(),
        spice_constellation_min_sats=_priority_spice_constellation_min_sats(),
        adaptive_enabled=_outer_route_adaptive_enabled(),
        adaptive_min_samples=_outer_route_min_samples(),
        failure_penalty_s=_outer_route_failure_penalty_s(),
        mc_process_min_samples=_outer_route_mc_process_min_samples(),
        mc_process_min_mission_s=_outer_route_mc_process_min_mission_s(),
        trace=_outer_route_trace_enabled()
    )
end

@inline function _outer_route_features(
    case::BenchmarkCase;
    spec::Union{Nothing, ProfileSpec}=nothing
)::ParallelProfiles.OuterRouteFeatures
    profile_name = isnothing(spec) ? lowercase(strip(get(ENV, "SPACEAGORA_PERF_PROFILE", "quick"))) : spec.name
    n_sats = length(case.args_template.dynamics_model.spacecraft)
    n_links = _total_links(case)
    max_links_per_sat = _max_links_per_sat(case)
    mission_time_s = case.args_template.mission_configuration.mission_time
    dt_max_orbit_s = case.args_template.integration_tolerances.dt_max_orbit
    dynamic_effectors = case.args_template.dynamics_model.dynamic_effectors
    density_model = case.args_template.environment_model.density_model
    has_nbody = _has_nbody_dynamic_effector(dynamic_effectors)
    has_srp = _has_srp_dynamic_effector(dynamic_effectors)
    harmonics_degree = _max_harmonics_degree(dynamic_effectors)
    has_control = !isempty(case.args_template.control_model.control_effectors)
    control_effector_count = length(case.args_template.control_model.control_effectors)
    orientation_on = case.args_template.mission_configuration.orientation_sim
    density_family = _density_model_family(density_model)
    solver_mode = _solver_mode_for_outer_route(case, profile_name)
    control_rate_s = _min_positive_rate(case.args_template.control_model.control_rates)
    guidance_rate_s = _min_positive_rate(case.args_template.guidance_model.guidance_rates)
    navigation_rate_s = _min_positive_rate(case.args_template.navigation_model.navigation_rates)
    gram_surrogate_enabled = _gram_surrogate_flag(density_model)
    gram_static_grid_enabled = _gram_static_grid_flag_for_case(case, density_model)
    thermal_enabled = _has_atmo_dynamic_effector(dynamic_effectors) || !(density_model isa NoAtmosphereModel)
    dynamic_effector_count = length(dynamic_effectors)
    effector_cost_class = _dynamic_effector_cost_class(dynamic_effectors, control_effector_count)
    mc_samples = isnothing(spec) ? 0 : spec.montecarlo_samples
    return ParallelProfiles.OuterRouteFeatures(
        category=case.category,
        n_sats=n_sats,
        n_links=n_links,
        max_links_per_sat=max_links_per_sat,
        mission_time_s=mission_time_s,
        has_nbody=has_nbody,
        has_srp=has_srp,
        harmonics_degree=harmonics_degree,
        has_control=has_control,
        orientation_on=orientation_on,
        density_family=density_family,
        solver_mode=solver_mode,
        dt_max_orbit_s=dt_max_orbit_s,
        control_rate_s=control_rate_s,
        guidance_rate_s=guidance_rate_s,
        navigation_rate_s=navigation_rate_s,
        gram_surrogate_enabled=gram_surrogate_enabled,
        gram_static_grid_enabled=gram_static_grid_enabled,
        control_effector_count=control_effector_count,
        thermal_enabled=thermal_enabled,
        dynamic_effector_count=dynamic_effector_count,
        effector_cost_class=effector_cost_class,
        montecarlo_samples=mc_samples
    )
end

@inline function _priority_outer_route_montecarlo(case::BenchmarkCase, spec::Union{Nothing, ProfileSpec}=nothing)::Symbol
    features = _outer_route_features(case; spec=spec)
    return ParallelProfiles.default_outer_route(
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function _outer_route_candidates(case::BenchmarkCase; spec::Union{Nothing, ProfileSpec}=nothing)::Vector{Symbol}
    features = _outer_route_features(case; spec=spec)
    return ParallelProfiles.outer_route_candidates(
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function _outer_route_signature(case::BenchmarkCase)::String
    return ParallelProfiles.outer_route_signature(_outer_route_features(case))
end

@inline function _outer_route_stats_snapshot(signature::String)
    return ParallelProfiles.outer_route_stats_snapshot(_PERF_OUTER_ROUTE_STATE[], signature)
end

function _record_outer_route_feedback!(
    case::BenchmarkCase,
    rows::AbstractVector{<:NamedTuple};
    route::Symbol
)
    features = _outer_route_features(case)
    successes = 0
    failures = 0
    elapsed_success_s = 0.0
    for row in rows
        if !hasproperty(row, :solve_success)
            continue
        end
        if row.solve_success === true
            if hasproperty(row, :total_time_s) && row.total_time_s isa Real && isfinite(Float64(row.total_time_s))
                successes += 1
                elapsed_success_s += Float64(row.total_time_s)
            end
        elseif row.solve_success === false
            failures += 1
        end
    end
    ParallelProfiles.record_outer_route_feedback!(
        _PERF_OUTER_ROUTE_STATE[],
        features;
        route=route,
        successes=successes,
        failures=failures,
        elapsed_success_s=elapsed_success_s,
        tuning=_outer_route_tuning()
    )
    return nothing
end

@inline function _priority_outer_route(case::BenchmarkCase)::Symbol
    features = _outer_route_features(case)
    return ParallelProfiles.default_outer_route(
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function parallel_priority_plan(case::BenchmarkCase, outer_route::Symbol)::ParallelPriorityPlan
    resolved_outer = outer_route == :auto ? _priority_outer_route(case) : outer_route
    n_sats = length(case.args_template.dynamics_model.spacecraft)
    n_links = _total_links(case)
    dynamic_effectors = case.args_template.dynamics_model.dynamic_effectors
    has_atmo = _has_atmo_dynamic_effector(dynamic_effectors)
    has_nbody = _has_nbody_dynamic_effector(dynamic_effectors)
    has_density = has_atmo || !(case.args_template.environment_model.density_model isa NoAtmosphereModel)
    has_control = !isempty(case.args_template.control_model.control_effectors)
    sat_threshold = _priority_inner_sat_threshold()
    link_threshold = _priority_inner_link_threshold()

    if resolved_outer != :none
        return ParallelPriorityPlan(
            outer_route=resolved_outer,
            density_mode="off",
            control_mode="off",
            multibody_mode="off"
        )
    end

    density_mode = has_density && n_sats >= sat_threshold ? "auto" : "off"
    control_mode = has_control && n_sats >= sat_threshold ? "auto" : "off"
    multibody_mode = (has_atmo || has_nbody) && n_links >= link_threshold ? "auto" : "off"

    return ParallelPriorityPlan(
        outer_route=resolved_outer,
        density_mode=density_mode,
        control_mode=control_mode,
        multibody_mode=multibody_mode
    )
end

@inline function _existing_or_policy_env(name::String, value::String)::String
    # Only respect user-provided baseline overrides captured at startup.
    # This avoids inheriting transient values left by other benchmark phases.
    baseline = get(_PERF_POLICY_ENV_BASELINE, name, nothing)
    if baseline === nothing
        return value
    end
    token = lowercase(strip(baseline))
    valid = if name == "SPACEAGORA_OUTER_PARALLEL_ACTIVE"
        token in ("1", "0", "true", "false", "yes", "no", "on", "off")
    else
        token in ("off", "none", "serial", "0", "false", "no", "on", "thread", "threads", "1", "true", "yes", "auto")
    end
    return valid ? baseline : value
end

@inline function parallel_priority_env_pairs(plan::ParallelPriorityPlan)::Vector{Pair{String, String}}
    outer_flag = plan.outer_route == :none ? "0" : "1"
    return Pair{String, String}[
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _existing_or_policy_env("SPACEAGORA_OUTER_PARALLEL_ACTIVE", outer_flag),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL", plan.density_mode),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL", plan.control_mode),
        "SPACEAGORA_MULTIBODY_PARALLEL" => _existing_or_policy_env("SPACEAGORA_MULTIBODY_PARALLEL", plan.multibody_mode),
    ]
end

@inline function _merged_env_pairs(
    base_pairs::Vector{Pair{String, String}},
    override_pairs::Vector{Pair{String, String}}
)::Vector{Pair{String, String}}
    isempty(override_pairs) && return base_pairs
    merged = copy(base_pairs)
    for pair in override_pairs
        key = first(pair)
        idx = findfirst(p -> first(p) == key, merged)
        if idx === nothing
            push!(merged, key => last(pair))
        else
            merged[idx] = key => last(pair)
        end
    end
    return merged
end

@inline function case_env_pairs(
    case::BenchmarkCase,
    plan::ParallelPriorityPlan
)::Vector{Pair{String, String}}
    return _merged_env_pairs(parallel_priority_env_pairs(plan), case.env_overrides)
end

@inline function _env_pairs_key(env_pairs::Vector{Pair{String, String}})::String
    isempty(env_pairs) && return ""
    return join([string(first(p), "=", last(p)) for p in env_pairs], ";")
end

function _thread_plan_groups(
    cases::Vector{BenchmarkCase},
    indices::Vector{Int},
    outer_route::Symbol
)::Vector{Tuple{Vector{Pair{String, String}}, Vector{Tuple{Int, ParallelPriorityPlan}}}}
    grouped_pairs = Dict{String, Vector{Pair{String, String}}}()
    grouped_payload = Dict{String, Vector{Tuple{Int, ParallelPriorityPlan}}}()
    ordered_keys = String[]
    for idx in indices
        plan = parallel_priority_plan(cases[idx], outer_route)
        env_pairs = case_env_pairs(cases[idx], plan)
        key = _env_pairs_key(env_pairs)
        if !haskey(grouped_payload, key)
            grouped_pairs[key] = env_pairs
            grouped_payload[key] = Tuple{Int, ParallelPriorityPlan}[]
            push!(ordered_keys, key)
        end
        push!(grouped_payload[key], (idx, plan))
    end
    return [(grouped_pairs[key], grouped_payload[key]) for key in ordered_keys]
end

function auto_backend_for_case(
    case::BenchmarkCase;
    spec::Union{Nothing, ProfileSpec}=nothing
)::Symbol
    features = _outer_route_features(case; spec=spec)
    return ParallelProfiles.select_outer_route!(
        _PERF_OUTER_ROUTE_STATE[],
        features;
        tuning=_outer_route_tuning(),
        machine_class=_machine_parallel_class(),
        threads_available=perf_threads_available(),
        parallel_enabled=perf_parallel_enabled()
    )
end

@inline function perf_process_workers_target()::Int
    raw = strip(get(ENV, "SPACEAGORA_PERF_PROCS", "0"))
    if isempty(raw) || raw == "0"
        return max(1, Sys.CPU_THREADS - 1)
    end
    n = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_PERF_PROCS must be an integer, got '$raw'"))
    end
    return max(1, n)
end

@inline function _is_valid_worker_project(path::String)::Bool
    return isdir(path) && (
        isfile(joinpath(path, "Project.toml")) ||
        isfile(joinpath(path, "JuliaProject.toml"))
    )
end

@inline function _project_candidate_path(raw::AbstractString)::String
    token = strip(String(raw))
    isempty(token) && return ""
    return normpath(isabspath(token) ? token : joinpath(REPO_ROOT, token))
end

function _resolve_perf_worker_project_path()::String
    checked = String[]
    override_raw = strip(get(ENV, "SPACEAGORA_PERF_WORKER_PROJECT", ""))
    if !isempty(override_raw)
        override_path = _project_candidate_path(override_raw)
        push!(checked, override_path)
        if _is_valid_worker_project(override_path)
            return override_path
        end
    end

    for candidate in (joinpath(REPO_ROOT, ".AGORA"), REPO_ROOT)
        candidate_path = normpath(candidate)
        candidate_path in checked || push!(checked, candidate_path)
        if _is_valid_worker_project(candidate_path)
            return candidate_path
        end
    end

    throw(ArgumentError(
        "Unable to resolve process-worker Julia project. " *
        "Set SPACEAGORA_PERF_WORKER_PROJECT to a valid project directory " *
        "(must contain Project.toml or JuliaProject.toml). " *
        "Checked: $(join(checked, ", "))."
    ))
end

const _perf_workers_initialized = Ref(false)
const _perf_worker_planet_cache = Ref{Any}(nothing)
const _perf_worker_mars_cache = Ref{Any}(nothing)
const _perf_harmonics20_cache = Ref{Any}(nothing)
const _perf_harmonics50_cache = Ref{Any}(nothing)
const _perf_mars_gram_point_density_cache = Ref{Any}(nothing)

function ensure_perf_workers!()
    _perf_workers_initialized[] && return nothing
    target_workers = perf_process_workers_target()
    missing_workers = target_workers - nworkers()
    if missing_workers > 0
        worker_project = _resolve_perf_worker_project_path()
        addprocs(
            missing_workers;
            exeflags=`--startup-file=no --project=$(worker_project)`
        )
    end
    script_path = abspath(@__FILE__)
    @everywhere workers() include($script_path)
    for w in workers()
        remotecall_wait(perf_worker_planet, w)
        remotecall_wait(perf_worker_mars, w)
        remotecall_wait(perf_worker_prime_gram_bindings!, w)
    end
    _perf_workers_initialized[] = true
    return nothing
end

function perf_worker_planet()::Earth
    cached = _perf_worker_planet_cache[]
    if cached === nothing
        cached = Earth("", SPICE_PATH)
        _perf_worker_planet_cache[] = cached
    end
    return cached::Earth
end

function perf_worker_mars()::Mars
    cached = _perf_worker_mars_cache[]
    if cached === nothing
        cached = Mars("", SPICE_PATH)
        _perf_worker_mars_cache[] = cached
    end
    return cached::Mars
end

function perf_worker_prime_gram_bindings!()::Nothing
    # Julia 1.12 is stricter about global bindings during cross-process
    # deserialization; eagerly load both GRAM wrappers on each worker.
    _build_earth_gram_point_density()
    _build_mars_gram_point_density()
    return nothing
end

function _perf_harmonics20_model(planet::Earth)
    cached = _perf_harmonics20_cache[]
    if cached === nothing
        cached = GravitationalHarmonicsModel(20, 20, EARTH_HARMONICS_FILE, planet)
        _perf_harmonics20_cache[] = cached
    end
    return deepcopy(cached)
end

function _perf_harmonics50_model(planet::Earth)
    cached = _perf_harmonics50_cache[]
    if cached === nothing
        cached = GravitationalHarmonicsModel(50, 50, EARTH_HARMONICS_FILE, planet)
        _perf_harmonics50_cache[] = cached
    end
    return deepcopy(cached)
end

function _perf_mars_gram_point_density_model()
    cached = _perf_mars_gram_point_density_cache[]
    if cached === nothing
        cached = _build_mars_gram_point_density()
        _perf_mars_gram_point_density_cache[] = cached
    end
    return deepcopy(cached)
end

@inline function orbital_period_seconds(spacecraft::SpacecraftModel, planet::AbstractPlanet)::Float64
    a = spacecraft.initial_condition.a
    if !isfinite(a) || a <= 0.0
        throw(ArgumentError("Invalid semimajor axis for orbital period calculation: $a"))
    end
    return 2π * sqrt(a^3 / planet.μ)
end

function _entry_orbital_elements_from_gamma_v_h(
    planet::AbstractPlanet;
    gamma_deg::Float64,
    v_mps::Float64,
    h_m::Float64
)::NamedTuple{(:a, :e, :ν_deg), Tuple{Float64, Float64, Float64}}
    v_mps > 0.0 || throw(ArgumentError("Entry speed must be > 0 m/s (got $v_mps)."))
    h_m >= 0.0 || throw(ArgumentError("Entry altitude h must be >= 0 m (got $h_m)."))

    r = planet.Rp_e + h_m
    γ = deg2rad(gamma_deg)
    μ = planet.μ
    specific_energy = 0.5 * v_mps^2 - μ / r
    abs(specific_energy) > 1e-12 || throw(ArgumentError("Parabolic entry state is unsupported (specific energy ~ 0)."))

    a = -μ / (2.0 * specific_energy)
    h_spec = r * v_mps * cos(γ)
    p = h_spec^2 / μ
    e_sq = 1.0 - p / a
    e_sq >= -1e-10 || throw(ArgumentError("Invalid entry state (computed e^2=$e_sq < 0). Check gamma/v/h."))
    e = sqrt(max(0.0, e_sq))

    ν_deg = if e <= 1e-10
        180.0
    else
        cos_ν = clamp((p / r - 1.0) / e, -1.0, 1.0)
        rdot = v_mps * sin(γ)
        sin_ν = clamp((rdot * h_spec) / (μ * e), -1.0, 1.0)
        rad2deg(mod(atan(sin_ν, cos_ν), 2π))
    end

    return (a=a, e=e, ν_deg=ν_deg)
end

function make_spacecraft(
    planet::AbstractPlanet;
    id::Int=1,
    ra_alt_m::Float64=550e3,
    rp_alt_m::Float64=500e3,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=170.0,
    with_panel::Bool=true,
    panel_count::Int=1,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing,
    root_mass::Float64=500.0,
    root_area::Float64=12.0,
    panel_mass::Float64=30.0,
    panel_area::Float64=6.0,
    panel_offset_y::Float64=1.2,
    prop_mass::Float64=0.0
)::SpacecraftModel
    root = Link{0}(root=true, m=root_mass, ref_area=root_area)
    links = Link[root]
    if with_panel
        panel_count >= 1 || throw(ArgumentError("panel_count must be >= 1 when with_panel=true (got $panel_count)."))
        if panel_count == 1
            panel = Link{0}(root=false, m=panel_mass, ref_area=panel_area, r=MVector{3, Float64}(0.0, panel_offset_y, 0.0))
            push!(links, panel)
        else
            # Spread appended bodies around the root to keep a balanced 5-body orientation benchmark.
            for panel_idx in 1:panel_count
                θ = 2π * (panel_idx - 1) / panel_count
                panel = Link{0}(
                    root=false,
                    m=panel_mass,
                    ref_area=panel_area,
                    r=MVector{3, Float64}(panel_offset_y * cos(θ), panel_offset_y * sin(θ), 0.0)
                )
                push!(links, panel)
            end
        end
    end

    ra = planet.Rp_e + ra_alt_m
    rp = planet.Rp_e + rp_alt_m

    ic = if isnothing(orientation_state)
        InitialCondition(ra=ra, rp=rp, i=i_deg, ω=ω_deg, Ω=Ω_deg, ν=ν_deg)
    else
        q0, w0 = orientation_state
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = sum(link.m for link in links)
    return SpacecraftModel(
        Joint[],
        links,
        root,
        true,
        dry_mass,
        prop_mass,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end

function make_blunted_cone_entry_spacecraft(
    planet::AbstractPlanet;
    id::Int=1,
    gamma_deg::Float64=-11.5,
    v_mps::Float64=5500.0,
    h_m::Float64=130e3,
    i_deg::Float64=51.6,
    ω_deg::Float64=30.0,
    Ω_deg::Float64=25.0,
    root_mass::Float64=320.0,
    base_radius_m::Float64=0.89,
    body_length_m::Float64=1.2,
    reflection_coefficient::Float64=1.0,
    prop_mass::Float64=0.0
)::SpacecraftModel
    base_radius_m > 0.0 || throw(ArgumentError("base_radius_m must be > 0 (got $base_radius_m)."))
    body_length_m > 0.0 || throw(ArgumentError("body_length_m must be > 0 (got $body_length_m)."))

    entry_oe = _entry_orbital_elements_from_gamma_v_h(
        planet;
        gamma_deg=gamma_deg,
        v_mps=v_mps,
        h_m=h_m
    )

    root = Link{0}(
        root=true,
        m=root_mass,
        ref_area=π * base_radius_m^2,
        dims=MVector{3, Float64}(body_length_m, 2.0 * base_radius_m, 2.0 * base_radius_m),
        reflection_coefficient=reflection_coefficient
    )
    links = Link[root]
    ic = InitialCondition(entry_oe.a, entry_oe.e, i_deg, ω_deg, Ω_deg, entry_oe.ν_deg)
    dry_mass = sum(link.m for link in links)

    return SpacecraftModel(
        Joint[],
        links,
        root,
        true,
        dry_mass,
        prop_mass,
        root.inertia,
        0,
        0,
        ic,
        id
    )
end

function make_constellation(
    planet::AbstractPlanet,
    n::Int;
    with_panel::Bool=false,
    panel_count::Int=1
)::Vector{SpacecraftModel}
    sats = SpacecraftModel[]
    for i in 1:n
        ra_alt = 540e3 + 20e3 * (i - 1)
        rp_alt = 470e3 + 15e3 * (i - 1)
        if rp_alt >= ra_alt
            ra_alt = rp_alt + 8e3
        end
        ν = 140.0 + 180.0 * (i - 1) / n
        push!(
            sats,
            make_spacecraft(
                planet;
                id=i,
                ra_alt_m=ra_alt,
                rp_alt_m=rp_alt,
                ν_deg=ν,
                with_panel=with_panel,
                panel_count=panel_count
            )
        )
    end
    return sats
end

function build_config(;
    planet::AbstractPlanet,
    spacecraft::Vector{SpacecraftModel},
    mission_time_s::Float64,
    orientation_sim::Bool,
    dynamic_effectors::Tuple,
    mission_type::MissionType=MissionTime,
    mission_keplerian::Bool=true,
    mission_orbits::Int=1,
    density_model=NoAtmosphereModel(),
    guidance_effectors::Tuple=(),
    guidance_rates::Vector{Float64}=Float64[],
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    dt_max_orbit::Float64=1.0,
    reltol_orbit::Float64=1e-9,
    abstol_orbit::Float64=1e-9
)::SimulationConfiguration
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false,
            verbose=false,
            generate_plots=false,
            normalize=false,
            save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=mission_type,
            keplerian=mission_keplerian,
            number_of_orbits=mission_orbits,
            mission_time=mission_time_s,
            orientation_sim=orientation_sim,
            num_steps_to_save=400
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=120.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=guidance_effectors, guidance_rates=guidance_rates),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=reltol_orbit,
            abstol_orbit=abstol_orbit,
            dt_max_orbit=dt_max_orbit
        )
    )
end

Base.@kwdef struct _GRAMOfflineSurrogateFileBase
    planet_name::String = "earth"
end

function SimulationModel.EnvironmentModels._gram_point_density(
    model::_GRAMOfflineSurrogateFileBase,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    # Offline benchmark path: keep density lookup file-backed and avoid native point GRAM calls.
    return 0.0, 200.0, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function _build_earth_gram_surrogate_density()
    isfile(EARTH_GRAM_SURROGATE_FILE) || throw(ArgumentError("GRAM surrogate file not found: $(EARTH_GRAM_SURROGATE_FILE)"))
    base_model = _GRAMOfflineSurrogateFileBase(planet_name="earth")
    return GRAMAtmosphereModelSurrogate(base_model, EARTH_GRAM_SURROGATE_FILE, nothing)
end

function _build_earth_gram_point_density()
    gram_root = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0")
    # GRAM Python bindings can be sensitive to world-age in long-lived Julia sessions.
    return Base.invokelatest(
        GRAMAtmosphereModel;
        gram_directory=gram_root,
        gram_data_directory=gram_root,
        spice_directory=SPICE_PATH,
        planet_name="earth",
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
    )
end

function _build_mars_gram_point_density()
    gram_root = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0")
    # GRAM Python bindings can be sensitive to world-age in long-lived Julia sessions.
    return Base.invokelatest(
        GRAMAtmosphereModel;
        gram_directory=gram_root,
        gram_data_directory=gram_root,
        spice_directory=SPICE_PATH,
        planet_name="mars",
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0)
    )
end

function make_montecarlo_config(seed::Int, planet::Earth, mission_time_s::Float64)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 530e3 + randn(rng) * 20e3
    rp_alt = 490e3 + randn(rng) * 20e3
    rp_alt = max(rp_alt, 120e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 8e3
    end

    sc = make_spacecraft(
        planet;
        id=1,
        with_panel=false,
        root_mass=160.0,
        root_area=1.5,
        prop_mass=20.0,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=28.0 + randn(rng) * 0.8,
        ω_deg=15.0 + randn(rng) * 2.0,
        Ω_deg=20.0 + randn(rng) * 2.0,
        ν_deg=160.0 + randn(rng) * 8.0
    )

    thruster = BaseThrusterModel(
        thrust=[0.6 + rand(rng) * 0.5],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[0.0],
        stop_burn_time=[180.0 + rand(rng) * 40.0],
        Isp=[280.0 + rand(rng) * 40.0]
    )

    return build_config(
        planet=planet,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster,),
        control_rates=[1.0],
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

@inline function _montecarlo_scenario_catalog()::Vector{NamedTuple}
    return NamedTuple[
        (
            variant=:mars_aerobraking,
            name="montecarlo_mars_aerobraking",
            description="Mars aerobraking mission randomized per seed (aero + MarsGRAM point-to-point)"
        ),
        (
            variant=:multi_sat,
            name="montecarlo_multi_sat",
            description="4-spacecraft constellation randomized per seed"
        ),
        (
            variant=:high_accuracy,
            name="montecarlo_high_accuracy",
            description="High-accuracy single-spacecraft mission randomized per seed (L50 harmonics)"
        )
    ]
end

@inline function _active_montecarlo_scenarios()::Vector{NamedTuple}
    return _montecarlo_scenario_catalog()
end

@inline function _montecarlo_batch_mission_time_s(spec::ProfileSpec, variant::Symbol)::Float64
    if variant == :mars_aerobraking
        # Keep Mars aerobraking seeds long enough to include drag-passage behavior.
        return max(spec.montecarlo_mission_s, spec.mission_long_s)
    end
    return spec.montecarlo_mission_s
end

function make_montecarlo_mars_aerobraking_config(
    seed::Int,
    mars::Mars,
    mission_time_s::Float64
)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 4500e3 + randn(rng) * 220e3
    rp_alt = clamp(120e3 + randn(rng) * 18e3, 95e3, 180e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 120e3
    end

    sc = make_spacecraft(
        mars;
        id=1,
        with_panel=false,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=93.0 + randn(rng) * 0.8,
        ω_deg=80.0 + randn(rng) * 3.0,
        Ω_deg=30.0 + randn(rng) * 3.0,
        ν_deg=180.0 + randn(rng) * 10.0
    )

    return build_config(
        planet=mars,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        mission_keplerian=false,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        density_model=_perf_mars_gram_point_density_model(),
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

function make_montecarlo_multi_sat_config(
    seed::Int,
    planet::Earth,
    mission_time_s::Float64
)::SimulationConfiguration
    rng = MersenneTwister(seed)
    spacecraft = SpacecraftModel[]
    for i in 1:4
        ra_alt = 540e3 + 20e3 * (i - 1) + randn(rng) * 8e3
        rp_alt = 470e3 + 15e3 * (i - 1) + randn(rng) * 8e3
        rp_alt = max(rp_alt, 120e3)
        if rp_alt >= ra_alt
            ra_alt = rp_alt + 8e3
        end
        ν = 140.0 + 180.0 * (i - 1) / 4 + randn(rng) * 5.0
        push!(
            spacecraft,
            make_spacecraft(
                planet;
                id=i,
                with_panel=false,
                ra_alt_m=ra_alt,
                rp_alt_m=rp_alt,
                i_deg=35.0 + randn(rng) * 0.4,
                ω_deg=40.0 + randn(rng) * 1.0,
                Ω_deg=10.0 + randn(rng) * 1.0,
                ν_deg=ν
            )
        )
    end
    harmonics20 = _perf_harmonics20_model(planet)
    return build_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(), harmonics20),
        dt_max_orbit=1.0,
        reltol_orbit=1e-9,
        abstol_orbit=1e-9
    )
end

function make_montecarlo_high_accuracy_config(
    seed::Int,
    planet::Earth,
    mission_time_s::Float64
)::SimulationConfiguration
    rng = MersenneTwister(seed)
    ra_alt = 550e3 + randn(rng) * 12e3
    rp_alt = 500e3 + randn(rng) * 12e3
    rp_alt = max(rp_alt, 140e3)
    if rp_alt >= ra_alt
        ra_alt = rp_alt + 8e3
    end

    sc = make_spacecraft(
        planet;
        id=1,
        with_panel=false,
        ra_alt_m=ra_alt,
        rp_alt_m=rp_alt,
        i_deg=35.0 + randn(rng) * 0.3,
        ω_deg=40.0 + randn(rng) * 0.8,
        Ω_deg=10.0 + randn(rng) * 0.8,
        ν_deg=170.0 + randn(rng) * 3.0
    )
    harmonics50 = _perf_harmonics50_model(planet)
    return build_config(
        planet=planet,
        spacecraft=[sc],
        mission_time_s=mission_time_s,
        orientation_sim=false,
        dynamic_effectors=(InverseSquaredGravityModel(), harmonics50),
        dt_max_orbit=0.5,
        reltol_orbit=1e-10,
        abstol_orbit=1e-10
    )
end

function make_montecarlo_case(
    seed::Int,
    mission_time_s::Float64,
    variant::Symbol,
    planet::Earth;
    mars::Union{Nothing, Mars}=nothing
)::BenchmarkCase
    catalog = _active_montecarlo_scenarios()
    scenario_idx = findfirst(s -> s.variant == variant, catalog)
    scenario_idx === nothing && throw(ArgumentError("Unsupported Monte Carlo scenario variant '$variant'."))
    scenario_meta = catalog[scenario_idx]

    args_template = if variant == :mars_aerobraking
        mars_planet = isnothing(mars) ? perf_worker_mars() : mars
        make_montecarlo_mars_aerobraking_config(seed, mars_planet, mission_time_s)
    elseif variant == :multi_sat
        make_montecarlo_multi_sat_config(seed, planet, mission_time_s)
    elseif variant == :high_accuracy
        make_montecarlo_high_accuracy_config(seed, planet, mission_time_s)
    else
        throw(ArgumentError("Unsupported Monte Carlo scenario variant '$variant'."))
    end

    return BenchmarkCase(
        name=scenario_meta.name,
        category="montecarlo",
        description=scenario_meta.description,
        args_template=args_template,
        run_in_quick=true
    )
end

function build_cases(spec::ProfileSpec, planet::Earth)::Vector{BenchmarkCase}
    harmonics20 = GravitationalHarmonicsModel(20, 20, EARTH_HARMONICS_FILE, planet)
    harmonics50 = GravitationalHarmonicsModel(50, 50, EARTH_HARMONICS_FILE, planet)
    nbody_sun_moon = NBodyGravityModel(body_names=("Sun", "Moon"), primary_body_name="Earth", planet=planet)
    mars = Mars("", SPICE_PATH)

    q0 = normalize(SVector{4, Float64}(0.15, -0.05, 0.2, 0.96))
    w0 = SVector{3, Float64}(0.01, -0.02, 0.015)
    q1 = normalize(SVector{4, Float64}(0.12, -0.08, 0.24, 0.96))
    w1 = SVector{3, Float64}(0.012, -0.018, 0.02)

    sc_baseline = [make_spacecraft(planet; id=1, with_panel=false)]
    sc_orientation = [make_spacecraft(planet; id=1, with_panel=true, panel_count=4, orientation_state=(q0, w0))]
    sc_entry_shallow = [make_blunted_cone_entry_spacecraft(
        planet;
        id=1,
        gamma_deg=-8.5,
        v_mps=5200.0,
        h_m=135e3,
        root_mass=320.0,
        base_radius_m=0.89,
        body_length_m=1.2,
        prop_mass=0.0,
        i_deg=51.6,
        ω_deg=30.0,
        Ω_deg=25.0
    )]
    sc_entry_nominal = [make_blunted_cone_entry_spacecraft(
        planet;
        id=1,
        gamma_deg=-11.5,
        v_mps=5500.0,
        h_m=130e3,
        root_mass=320.0,
        base_radius_m=0.89,
        body_length_m=1.2,
        prop_mass=0.0,
        i_deg=51.6,
        ω_deg=30.0,
        Ω_deg=25.0
    )]
    sc_entry_steep = [make_blunted_cone_entry_spacecraft(
        planet;
        id=1,
        gamma_deg=-14.5,
        v_mps=5800.0,
        h_m=125e3,
        root_mass=320.0,
        base_radius_m=0.89,
        body_length_m=1.2,
        prop_mass=0.0,
        i_deg=51.6,
        ω_deg=30.0,
        Ω_deg=25.0
    )]
    sc_mars_aerobrake = [make_spacecraft(
        mars;
        id=1,
        with_panel=false,
        ra_alt_m=4500e3,
        rp_alt_m=120e3,
        i_deg=93.0,
        ω_deg=80.0,
        Ω_deg=30.0,
        ν_deg=180.0
    )]
    earth_gram_point_density = _build_earth_gram_point_density()
    mars_gram_point_density = _build_mars_gram_point_density()
    earth_gram_surrogate_density = _build_earth_gram_surrogate_density()
    thermal_stress_density = SimulationModel.EnvironmentModels.PolynomialFitAtmosphereModel([-27.0])
    multi_scaling_effectors = (InverseSquaredGravityModel(), harmonics20)
    sc_thermal_stress = make_constellation(planet, 8; with_panel=true, panel_count=12)
    sc_effector_stress6 = make_constellation(planet, 6; with_panel=false)
    sc_effector_stress12 = make_constellation(planet, 12; with_panel=false)
    effector_stress_effectors = (
        InverseSquaredGravityModel(),
        InverseSquaredJ2GravityModel(),
        SolarRadiationPressureModel(1.2, 16.0),
        SolarRadiationPressureModel(1.6, 24.0)
    )
    sc_proximity_fullstack = [
        make_spacecraft(
            planet;
            id=1,
            with_panel=true,
            orientation_state=(q0, w0),
            ra_alt_m=520e3,
            rp_alt_m=515e3,
            ν_deg=168.0,
            root_mass=280.0,
            root_area=4.0,
            panel_mass=12.0,
            panel_area=1.8,
            panel_offset_y=1.0,
            prop_mass=18.0
        ),
        make_spacecraft(
            planet;
            id=2,
            with_panel=true,
            orientation_state=(q1, w1),
            ra_alt_m=520.15e3,
            rp_alt_m=515.1e3,
            ν_deg=168.18,
            root_mass=280.0,
            root_area=4.0,
            panel_mass=12.0,
            panel_area=1.8,
            panel_offset_y=1.0,
            prop_mass=18.0
        )
    ]
    proximity_thruster = BaseThrusterModel(
        thrust=fill(0.28, 2),
        direction=fill(0.0, 2),
        Δv=fill(2.5, 2),
        start_burn_time=fill(-1.0, 2),
        stop_burn_time=fill(-1.0, 2),
        Isp=fill(290.0, 2)
    )
    proximity_guidance = (
        AerobrakingCampaignPropulsiveManeuverGuidanceModel(
            maneuver_orbit_number=[1],
            maneuver_Δv=[2.5]
        ),
        AerobrakingCampaignPropulsiveManeuverGuidanceModel(
            maneuver_orbit_number=[1],
            maneuver_Δv=[2.5]
        )
    )

    cases = BenchmarkCase[
        BenchmarkCase(
            name="single_orientation_aero",
            category="orientation",
            description="1 spacecraft (5-body), orientation dynamics on, aerodynamic model active",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_orientation,
                mission_time_s=spec.mission_short_s,
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM())
            )
        ),
        BenchmarkCase(
            name="thermal_8sat_panel12_aero",
            category="thermal_stress",
            description="8 spacecraft (13 links each) with aerodynamic model and fixed polynomial atmosphere to stress thermal callback throughput",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_thermal_stress,
                mission_time_s=min(spec.mission_short_s, 3600.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(thermal_stress_density),
                dt_max_orbit=0.5
            )
        ),
        BenchmarkCase(
            name="effector_6sat_dual_srp_stack",
            category="effector_stress",
            description="6 spacecraft with inverse gravity, J2, and dual SRP effectors (no atmosphere/control) to stress dynamic effector reduction with outer routing off (r2)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_effector_stress6,
                mission_time_s=min(spec.mission_short_s, 3600.0),
                orientation_sim=false,
                dynamic_effectors=effector_stress_effectors,
                density_model=NoAtmosphereModel(),
                dt_max_orbit=0.5
            )
        ),
        BenchmarkCase(
            name="effector_12sat_dual_srp_stack",
            category="effector_stress",
            description="12 spacecraft with inverse gravity, J2, and dual SRP effectors (no atmosphere/control) for larger multi-satellite effector scaling",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_effector_stress12,
                mission_time_s=min(spec.mission_short_s, 3600.0),
                orientation_sim=false,
                dynamic_effectors=effector_stress_effectors,
                density_model=NoAtmosphereModel(),
                dt_max_orbit=0.5
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_entry_earth_shallow",
            category="entry",
            description="1 blunted-cone entry spacecraft, Earth shallow entry from gamma/v/h (target 1 entry interface downcrossing, drag + aero, Earth GRAM point-to-point)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_entry_shallow,
                mission_time_s=max(spec.mission_short_s, 2400.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(earth_gram_point_density),
                dt_max_orbit=0.5
            ),
            entry_target_count_override=1
        ),
        BenchmarkCase(
            name="single_entry_earth_nominal",
            category="entry",
            description="1 blunted-cone entry spacecraft, Earth nominal entry from gamma/v/h (target 1 entry interface downcrossing, drag + aero, Earth GRAM point-to-point)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_entry_nominal,
                mission_time_s=max(spec.mission_short_s, 1800.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(earth_gram_point_density),
                dt_max_orbit=0.5
            ),
            entry_target_count_override=1
        ),
        BenchmarkCase(
            name="single_entry_earth_steep",
            category="entry",
            description="1 blunted-cone entry spacecraft, Earth steep entry from gamma/v/h (target 1 entry interface downcrossing, drag + aero, Earth GRAM point-to-point)",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_entry_steep,
                mission_time_s=max(spec.mission_short_s, 1200.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(earth_gram_point_density),
                dt_max_orbit=0.25
            ),
            entry_target_count_override=1
        ),
        BenchmarkCase(
            name="multi_4_gravity",
            category="satellite_scaling",
            description="4 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 4; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            )
        ),
        BenchmarkCase(
            name="multi_8_gravity",
            category="satellite_scaling",
            description="8 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 8; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_8_gravity_surrogate_cached",
            category="satellite_scaling",
            description="8 spacecraft, L20 harmonics with GRAM surrogate density and track-cache enabled",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 8; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false,
            env_overrides=Pair{String, String}[
                "SPACEAGORA_GRAM_TRACK_CACHE" => "on"
            ]
        ),
        BenchmarkCase(
            name="multi_16_gravity",
            category="satellite_scaling",
            description="16 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 16; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_32_gravity",
            category="satellite_scaling",
            description="32 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 32; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_64_gravity",
            category="satellite_scaling",
            description="64 spacecraft, L20 harmonics with GRAM surrogate density from file",
            args_template=build_config(
                planet=planet,
                spacecraft=make_constellation(planet, 64; with_panel=false),
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=multi_scaling_effectors,
                density_model=deepcopy(earth_gram_surrogate_density)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="single_j2",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square + J2 gravity",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(),)
            )
        ),
        BenchmarkCase(
            name="single_nbody_sun_moon",
            category="dynamics_fidelity",
            description="1 spacecraft, J2 gravity + N-body Sun/Moon perturbations",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(), nbody_sun_moon),
                dt_max_orbit=10.0
            )
        ),
        BenchmarkCase(
            name="single_harmonics_l20",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=20",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20)
            )
        ),
        BenchmarkCase(
            name="single_harmonics_l50",
            category="dynamics_fidelity",
            description="1 spacecraft, inverse-square gravity + spherical harmonics L=M=50",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_baseline,
                mission_time_s=spec.mission_short_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics50)
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="proximity_2sat_orientation_fullstack_gnc_highrate",
            category="rpo_gnc",
            description="2-spacecraft close-proximity operations with orientation on, high-rate guidance, and BaseThrusterModel control callback",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_proximity_fullstack,
                mission_time_s=min(spec.mission_short_s, 2400.0),
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(),),
                guidance_effectors=proximity_guidance,
                guidance_rates=[0.1, 0.1],
                control_effectors=(proximity_thruster,),
                control_rates=[0.1],
                dt_max_orbit=0.2
            )
        ),
        BenchmarkCase(
            name="single_baseline_long_mission",
            category="mission_length",
            description="1 spacecraft Mars aerobraking-style long mission (10 orbits, atmospheric drag + aero, MarsGRAM point-to-point)",
            args_template=build_config(
                planet=mars,
                spacecraft=sc_mars_aerobrake,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                mission_type=MissionOrbits,
                mission_keplerian=false,
                mission_orbits=10,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(mars_gram_point_density),
                dt_max_orbit=1.0
            )
        )
    ]

    return cases
end

@inline function _split_variant_case_name(base::AbstractString, solver::AbstractString)::String
    return string(base, "__split_imex_", lowercase(String(solver)))
end

@inline function _split_base_scenario_name(name::AbstractString)::String
    token = "__split_imex_"
    idx = findfirst(token, String(name))
    if idx === nothing
        return String(name)
    end
    start_idx = first(idx)
    return String(name)[1:(start_idx - 1)]
end

@inline function _split_is_variant_case(name::AbstractString)::Bool
    return occursin("__split_imex_", String(name))
end

function _split_rollout_benchmark_cases(cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    if !_split_rollout_enabled()
        return cases
    end
    target_names = Set(_split_rollout_case_names())
    split_solvers = _split_rollout_solver_variants()
    expanded = copy(cases)
    for case in cases
        if !(case.name in target_names)
            continue
        end
        for split_solver in split_solvers
            push!(expanded, BenchmarkCase(
                name=_split_variant_case_name(case.name, split_solver),
                category=case.category,
                description=string(case.description, " [split_imex:", split_solver, "]"),
                args_template=case.args_template,
                run_in_quick=case.run_in_quick,
                solver_mode_override="split_imex",
                split_imex_solver_override=split_solver,
                entry_target_count_override=case.entry_target_count_override,
                env_overrides=copy(case.env_overrides)
            ))
        end
    end
    return expanded
end

@inline function _multirate_variant_case_name(base::AbstractString)::String
    return string(base, "__multirate")
end

function _multirate_rollout_benchmark_cases(cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    if !_multirate_rollout_enabled()
        return cases
    end
    target_names = Set(_multirate_rollout_case_names())
    expanded = copy(cases)
    for case in cases
        if !(case.name in target_names)
            continue
        end
        push!(expanded, BenchmarkCase(
            name=_multirate_variant_case_name(case.name),
            category=case.category,
            description=string(case.description, " [multirate]"),
            args_template=case.args_template,
            run_in_quick=case.run_in_quick,
            solver_mode_override="multirate",
            split_imex_solver_override=nothing,
            entry_target_count_override=case.entry_target_count_override,
            env_overrides=copy(case.env_overrides)
        ))
    end
    return expanded
end

function measure_case(
    case::BenchmarkCase,
    profile_name::String,
    repeat_idx::Int;
    seed::Union{Missing, Int}=missing,
    attempt::Int=1,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing
)
    timestamp_utc = string(now(UTC))
    args_meta = case.args_template
    n_sats = length(args_meta.dynamics_model.spacecraft)
    mission_time_s = args_meta.mission_configuration.mission_time
    density_model = args_meta.environment_model.density_model
    density_family = _density_model_family(density_model)
    gram_surrogate_enabled = _gram_surrogate_flag(density_model)
    gram_static_grid_enabled = (
        density_model isa GRAMAtmosphereModel || density_model isa GRAMAtmosphereModelSurrogate
    ) && (
        (_env_bool_token(get(ENV, "SPACEAGORA_GRAM_STATIC_GRID", "off")) === true) ||
        _gram_static_grid_flag_for_case(case, density_model)
    )
    gram_track_cache_mode = _gram_track_cache_mode_for_case(case)
    density_backend_bucket = _density_backend_bucket(
        density_family,
        gram_surrogate_enabled,
        gram_static_grid_enabled,
        gram_track_cache_mode
    )
    resolved_plan = isnothing(plan) ? ParallelPriorityPlan() : plan
    hardware = _runtime_hardware_snapshot()

    GC.gc()
    copy_timed = @timed deepcopy(case.args_template)
    args_run = copy_timed.value

    GC.gc()
    solve_timed = @timed begin
        try
            result = _run_perf_simulation(
                args_run;
                return_solution=true,
                return_solver_metadata=true,
                profile_name=profile_name,
                solver_mode_override=case.solver_mode_override,
                split_imex_solver_override=case.split_imex_solver_override,
                entry_target_count_override=case.entry_target_count_override
            )
            (ok=true, result=result, err=nothing, bt=nothing)
        catch err
            if err isa InterruptException
                rethrow()
            end
            if _perf_strict_errors()
                rethrow(err)
            end
            (ok=false, result=nothing, err=err, bt=catch_backtrace())
        end
    end

    copy_time_s = copy_timed.time
    solve_time_s = solve_timed.time
    total_time_s = copy_time_s + solve_time_s

    copy_bytes_mb = copy_timed.bytes / 1024^2
    solve_bytes_mb = solve_timed.bytes / 1024^2
    total_bytes_mb = copy_bytes_mb + solve_bytes_mb

    copy_alloc_calls = _alloc_calls(copy_timed.gcstats)
    solve_alloc_calls = _alloc_calls(solve_timed.gcstats)

    solver_mode = missing
    solver_sequence = missing
    solver_fallback_used = missing
    solver_fallback_count = missing
    solver_fallback_trigger = missing
    solve_retcode = missing
    solve_success = false
    saved_points = missing
    accepted_steps = missing
    rejected_steps = missing
    sim_seconds_per_wall_second = missing
    satellite_sim_seconds_per_wall_second = missing
    policy_decisions_total = missing
    policy_threads_enabled_total = missing
    policy_density_threads_enabled = missing
    policy_control_threads_enabled = missing
    policy_multibody_threads_enabled = missing
    policy_other_threads_enabled = missing
    nbody_spkpos_runtime_calls = missing
    nbody_spkpos_cache_build_calls = missing
    nbody_spkpos_total_calls = missing
    srp_spkpos_runtime_calls = missing
    srp_spkpos_cache_build_calls = missing
    srp_spkpos_total_calls = missing
    planet_pxform_runtime_calls = missing
    planet_pxform_cache_build_calls = missing
    planet_pxform_total_calls = missing
    final_primary_pos_norm_m = missing
    final_primary_vel_norm_mps = missing
    final_primary_mass_kg = missing

    solve_payload = solve_timed.value
    if solve_payload.ok
        solve_result = solve_payload.result
        sol = solve_result.solution
        solver_mode = solve_result.solver_mode
        solver_trace = solve_result.solver_trace
        solver_sequence = isempty(solver_trace) ? missing : join([meta.solver for meta in solver_trace], "->")
        solver_fallback_count = count(meta -> meta.fallback_used, solver_trace)
        solver_fallback_used = solver_fallback_count > 0
        fallback_triggers = [meta.trigger_retcode for meta in solver_trace if !(meta.trigger_retcode isa Missing)]
        solver_fallback_trigger = isempty(fallback_triggers) ? missing : _safe_unique_join(fallback_triggers; delimiter="|")
        solve_retcode = string(sol.retcode)
        solve_success = _solve_success_for_case(sol, case)
        saved_points = length(sol.t)
        accepted_steps = _destat_int(sol, :naccept)
        rejected_steps = _destat_int(sol, :nreject)
        sim_seconds_per_wall_second = _safe_div(mission_time_s, total_time_s)
        satellite_sim_seconds_per_wall_second = _safe_div(mission_time_s * n_sats, total_time_s)
        terminal = _primary_terminal_state_metrics(sol)
        final_primary_pos_norm_m = terminal.pos_norm_m
        final_primary_vel_norm_mps = terminal.vel_norm_mps
        final_primary_mass_kg = terminal.mass_kg
        if hasproperty(solve_result, :parallel_policy)
            snapshot = solve_result.parallel_policy
            if !(snapshot isa Nothing)
                policy_decisions_total = getproperty(snapshot, :decisions_total)
                policy_threads_enabled_total = getproperty(snapshot, :threads_enabled_total)
                policy_density_threads_enabled = getproperty(snapshot, :density_threads_enabled)
                policy_control_threads_enabled = getproperty(snapshot, :control_threads_enabled)
                policy_multibody_threads_enabled = getproperty(snapshot, :multibody_threads_enabled)
                policy_other_threads_enabled = getproperty(snapshot, :other_threads_enabled)
            end
        end
        if hasproperty(solve_result, :spice_counters)
            counters = solve_result.spice_counters
            if !(counters isa Nothing)
                nbody_spkpos_runtime_calls = getproperty(counters, :nbody_spkpos_runtime_calls)
                nbody_spkpos_cache_build_calls = getproperty(counters, :nbody_spkpos_cache_build_calls)
                nbody_spkpos_total_calls = getproperty(counters, :nbody_spkpos_total_calls)
                srp_spkpos_runtime_calls = getproperty(counters, :srp_spkpos_runtime_calls)
                srp_spkpos_cache_build_calls = getproperty(counters, :srp_spkpos_cache_build_calls)
                srp_spkpos_total_calls = getproperty(counters, :srp_spkpos_total_calls)
                planet_pxform_runtime_calls = getproperty(counters, :planet_pxform_runtime_calls)
                planet_pxform_cache_build_calls = getproperty(counters, :planet_pxform_cache_build_calls)
                planet_pxform_total_calls = getproperty(counters, :planet_pxform_total_calls)
            end
        end
    else
        solve_err = solve_payload.err
        solve_bt = solve_payload.bt
        solve_retcode = _solve_retcode_from_error(solve_err)
        if ismissing(solve_retcode)
            solve_retcode = "Exception"
            @warn "[perf] case=$(case.name) repeat=$(repeat_idx) attempt=$(attempt) threw $(typeof(solve_err)): $(_perf_error_text(solve_err)) @ $(_perf_stack_head(solve_bt))"
        end
    end

    return (
        profile=profile_name,
        hardware_class=hardware.hardware_class,
        machine_label=hardware.machine_label,
        host_name=hardware.host_name,
        cpu_name=hardware.cpu_name,
        cpu_threads=hardware.cpu_threads,
        julia_threads=hardware.julia_threads,
        os=hardware.os,
        arch=hardware.arch,
        category=case.category,
        scenario=case.name,
        description=case.description,
        density_family=density_family,
        gram_surrogate_enabled=gram_surrogate_enabled,
        gram_static_grid_enabled=gram_static_grid_enabled,
        gram_track_cache_mode=gram_track_cache_mode,
        density_backend_bucket=density_backend_bucket,
        seed=seed,
        repeat=repeat_idx,
        attempt=attempt,
        satellites=n_sats,
        orientation=args_meta.mission_configuration.orientation_sim,
        mission_time_s=mission_time_s,
        outer_route=string(resolved_plan.outer_route),
        density_parallel_mode=resolved_plan.density_mode,
        control_parallel_mode=resolved_plan.control_mode,
        multibody_parallel_mode=resolved_plan.multibody_mode,
        dt_max_orbit_s=args_meta.integration_tolerances.dt_max_orbit,
        dynamic_effectors=_effector_signature(args_meta.dynamics_model.dynamic_effectors),
        control_effectors=_effector_signature(args_meta.control_model.control_effectors),
        copy_time_s=copy_time_s,
        solve_time_s=solve_time_s,
        total_time_s=total_time_s,
        copy_compile_time_s=copy_timed.compile_time,
        solve_compile_time_s=solve_timed.compile_time,
        copy_gctime_s=copy_timed.gctime,
        solve_gctime_s=solve_timed.gctime,
        solver_mode=solver_mode,
        split_imex_solver_override=isnothing(case.split_imex_solver_override) ? missing : String(case.split_imex_solver_override),
        solver_sequence=solver_sequence,
        solver_fallback_used=solver_fallback_used,
        solver_fallback_count=solver_fallback_count,
        solver_fallback_trigger=solver_fallback_trigger,
        solve_retcode=solve_retcode,
        solve_success=solve_success,
        copy_bytes_mb=copy_bytes_mb,
        solve_bytes_mb=solve_bytes_mb,
        total_bytes_mb=total_bytes_mb,
        copy_alloc_calls=copy_alloc_calls,
        solve_alloc_calls=solve_alloc_calls,
        saved_points=saved_points,
        accepted_steps=accepted_steps,
        rejected_steps=rejected_steps,
        policy_decisions_total=policy_decisions_total,
        policy_threads_enabled_total=policy_threads_enabled_total,
        policy_density_threads_enabled=policy_density_threads_enabled,
        policy_control_threads_enabled=policy_control_threads_enabled,
        policy_multibody_threads_enabled=policy_multibody_threads_enabled,
        policy_other_threads_enabled=policy_other_threads_enabled,
        nbody_spkpos_runtime_calls=nbody_spkpos_runtime_calls,
        nbody_spkpos_cache_build_calls=nbody_spkpos_cache_build_calls,
        nbody_spkpos_total_calls=nbody_spkpos_total_calls,
        srp_spkpos_runtime_calls=srp_spkpos_runtime_calls,
        srp_spkpos_cache_build_calls=srp_spkpos_cache_build_calls,
        srp_spkpos_total_calls=srp_spkpos_total_calls,
        planet_pxform_runtime_calls=planet_pxform_runtime_calls,
        planet_pxform_cache_build_calls=planet_pxform_cache_build_calls,
        planet_pxform_total_calls=planet_pxform_total_calls,
        final_primary_pos_norm_m=final_primary_pos_norm_m,
        final_primary_vel_norm_mps=final_primary_vel_norm_mps,
        final_primary_mass_kg=final_primary_mass_kg,
        sim_seconds_per_wall_second=sim_seconds_per_wall_second,
        satellite_sim_seconds_per_wall_second=satellite_sim_seconds_per_wall_second,
        timestamp_utc=timestamp_utc
    )
end

function run_warmup(case::BenchmarkCase, warmup::Int, profile_name::String)
    log_warmup = _perf_warmup_logs()
    for i in 1:warmup
        args_run = deepcopy(case.args_template)
        if log_warmup
            println("  warmup $(i)/$(warmup): start")
            flush(stdout)
        end
        warmup_started_ns = time_ns()
        try
            _run_perf_simulation(
                args_run;
                return_solution=false,
                profile_name=profile_name,
                solver_mode_override=case.solver_mode_override,
                split_imex_solver_override=case.split_imex_solver_override,
                entry_target_count_override=case.entry_target_count_override
            )
            warmup_elapsed_s = (time_ns() - warmup_started_ns) / 1e9
            if log_warmup
                println("  warmup $(i)/$(warmup): done total=$(round(warmup_elapsed_s; digits=3)) s")
                flush(stdout)
            end
        catch err
            if err isa InterruptException
                rethrow()
            end
            if _perf_strict_errors()
                rethrow(err)
            end
            err_bt = catch_backtrace()
            solve_retcode = _solve_retcode_from_error(err)
            warmup_elapsed_s = (time_ns() - warmup_started_ns) / 1e9
            if ismissing(solve_retcode)
                @warn "[warmup] $(case.name) $(i)/$(warmup) threw $(typeof(err)): $(_perf_error_text(err)) @ $(_perf_stack_head(err_bt)); continuing (elapsed=$(round(warmup_elapsed_s; digits=3)) s)"
            elseif log_warmup
                println("  warmup $(i)/$(warmup): failed retcode=$(solve_retcode), continuing (elapsed=$(round(warmup_elapsed_s; digits=3)) s)")
                flush(stdout)
            end
        end
    end
    return nothing
end

@inline function _case_sample_schedule(case::BenchmarkCase, spec::ProfileSpec)::Tuple{Int, Int}
    warmup = spec.warmup
    repeats = spec.repeats
    if case.category == "control_stress"
        if spec.name == "full"
            warmup = _control_stress_warmup_full()
            repeats = _control_stress_repeats_full()
        else
            warmup = min(warmup, 1)
            repeats = min(repeats, 2)
        end
    end
    return warmup, repeats
end

function _run_case_batch_core!(
    rows::Vector{NamedTuple},
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int,
    plan::ParallelPriorityPlan
)
    warmup_count, repeat_count = _case_sample_schedule(case, spec)
    println(
        "[$idx/$total] $(case.name) :: warmup x$(warmup_count), repeats x$(repeat_count), " *
        "outer=$(plan.outer_route), density=$(plan.density_mode), control=$(plan.control_mode), multibody=$(plan.multibody_mode), " *
        "solver_override=$(isnothing(case.solver_mode_override) ? "none" : case.solver_mode_override), " *
        "split_override=$(isnothing(case.split_imex_solver_override) ? "none" : case.split_imex_solver_override)"
    )
    run_warmup(case, warmup_count, spec.name)
    for rep in 1:repeat_count
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, rep; attempt=attempt, plan=plan)
            last_row = row
            if row.solve_success
                push!(rows, row)
                println("  repeat $(rep)/$(repeat_count) attempt $(attempt)/$(spec.max_attempts): total=$(round(row.total_time_s; digits=3)) s, solve=$(round(row.solve_time_s; digits=3)) s")
                flush(stdout)
                break
            end
            println("  repeat $(rep)/$(repeat_count) attempt $(attempt)/$(spec.max_attempts): failed retcode=$(row.solve_retcode), retrying")
            flush(stdout)
        end
        if !(last_row === nothing) && !last_row.solve_success
            push!(rows, last_row)
            println("  repeat $(rep)/$(repeat_count): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)")
            flush(stdout)
        end
    end
    return nothing
end

function run_case_batch!(
    rows::Vector{NamedTuple},
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int;
    outer_route::Symbol=:none,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing,
    apply_env::Bool=true
)
    resolved_plan = isnothing(plan) ? parallel_priority_plan(case, outer_route) : plan
    if apply_env
        env_pairs = case_env_pairs(case, resolved_plan)
        return withenv(env_pairs...) do
            _run_case_batch_core!(rows, case, spec, idx, total, resolved_plan)
        end
    end
    return _run_case_batch_core!(rows, case, spec, idx, total, resolved_plan)
end

function measure_montecarlo_seed(
    spec::ProfileSpec,
    planet::Earth,
    mission_time_s::Float64,
    seed::Int;
    variant::Symbol=:high_accuracy,
    mars::Union{Nothing, Mars}=nothing,
    outer_route::Symbol=:none,
    plan::Union{Nothing, ParallelPriorityPlan}=nothing,
    apply_env::Bool=true
)
    case = make_montecarlo_case(seed, mission_time_s, variant, planet; mars=mars)
    resolved_plan = isnothing(plan) ? parallel_priority_plan(case, outer_route) : plan
    run_seed = () -> begin
        last_row = nothing
        for attempt in 1:spec.max_attempts
            row = measure_case(case, spec.name, 1; seed=seed, attempt=attempt, plan=resolved_plan)
            last_row = row
            if row.solve_success
                return row, nothing
            end
        end
        return last_row, last_row === nothing ? "failed without attempt data" : "failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)"
    end
    if apply_env
        env_pairs = parallel_priority_env_pairs(resolved_plan)
        return withenv(env_pairs...) do
            run_seed()
        end
    end
    return run_seed()
end

function perf_worker_montecarlo_warmup(
    spec::ProfileSpec,
    mission_time_s::Float64,
    seed::Int,
    variant::Symbol,
    outer_route::Symbol=:process
)
    planet = perf_worker_planet()
    mars = perf_worker_mars()
    warmup_case = make_montecarlo_case(seed, mission_time_s, variant, planet; mars=mars)
    plan = parallel_priority_plan(warmup_case, outer_route)
    env_pairs = parallel_priority_env_pairs(plan)
    withenv(env_pairs...) do
        run_warmup(warmup_case, spec.warmup, spec.name)
    end
    return nothing
end

function perf_worker_measure_montecarlo_seed(
    spec::ProfileSpec,
    mission_time_s::Float64,
    seed::Int,
    variant::Symbol,
    outer_route::Symbol=:process
)
    planet = perf_worker_planet()
    mars = perf_worker_mars()
    return measure_montecarlo_seed(
        spec,
        planet,
        mission_time_s,
        seed;
        variant=variant,
        mars=mars,
        outer_route=outer_route
    )
end

function perf_worker_run_case_batch(
    case::BenchmarkCase,
    spec::ProfileSpec,
    idx::Int,
    total::Int,
    outer_route::Symbol=:process
)
    local_rows = NamedTuple[]
    run_case_batch!(local_rows, case, spec, idx, total; outer_route=outer_route)
    return local_rows
end

function perf_worker_measure_per_orbit_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    period_s::Float64,
    orbit_counts::Vector{Int},
    outer_route::Symbol=:process
)
    return measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts; outer_route=outer_route)
end

function run_montecarlo_batch!(rows::Vector{NamedTuple}, spec::ProfileSpec, planet::Earth)
    seeds = collect(1001:(1000 + spec.montecarlo_samples))
    scenarios = _active_montecarlo_scenarios()
    scenario_names = join([s.name for s in scenarios], ", ")
    println("[montecarlo] warmup x$(spec.warmup), seeds=$(length(seeds)), scenarios=$(scenario_names)")

    backend = perf_parallel_backend()
    mars = perf_worker_mars()

    for scenario in scenarios
        variant = scenario.variant
        mission_time_s = _montecarlo_batch_mission_time_s(spec, variant)
        warmup_case = make_montecarlo_case(first(seeds), mission_time_s, variant, planet; mars=mars)
        println("  scenario $(scenario.name) (mission_time=$(round(mission_time_s; digits=1)) s)")

        mc_backend = backend == :auto ? auto_backend_for_case(warmup_case; spec=spec) : backend
        if mc_backend == :threads && !_case_outer_threads_safe(warmup_case)
            mc_backend = :none
        end
        if mc_backend == :process
            ensure_perf_workers!()
            warmup_seed = first(seeds)
            for w in workers()
                remotecall_wait(perf_worker_montecarlo_warmup, w, spec, mission_time_s, warmup_seed, variant, :process)
            end
        else
            plan = parallel_priority_plan(warmup_case, mc_backend)
            env_pairs = parallel_priority_env_pairs(plan)
            withenv(env_pairs...) do
                run_warmup(warmup_case, spec.warmup, spec.name)
            end
        end

        seed_rows = Vector{NamedTuple}(undef, length(seeds))
        seed_msgs = Vector{String}(undef, length(seeds))

        if mc_backend == :process
            seed_results = pmap(seed -> perf_worker_measure_montecarlo_seed(spec, mission_time_s, seed, variant, :process), seeds)
            for i in eachindex(seeds)
                seed = seeds[i]
                row, err = seed_results[i]
                seed_rows[i] = row
                if row.solve_success
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
                else
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): $(err)"
                end
            end
        elseif mc_backend == :threads
            threaded_plan = parallel_priority_plan(warmup_case, :threads)
            threaded_env = parallel_priority_env_pairs(threaded_plan)
            withenv(threaded_env...) do
                Threads.@threads for i in eachindex(seeds)
                    seed = seeds[i]
                    row, err = measure_montecarlo_seed(
                        spec,
                        planet,
                        mission_time_s,
                        seed;
                        variant=variant,
                        mars=mars,
                        outer_route=:threads,
                        plan=threaded_plan,
                        apply_env=false
                    )
                    seed_rows[i] = row
                    if row.solve_success
                        seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
                    else
                        seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): $(err)"
                    end
                end
            end
        else
            for i in eachindex(seeds)
                seed = seeds[i]
                row, err = measure_montecarlo_seed(
                    spec,
                    planet,
                    mission_time_s,
                    seed;
                    variant=variant,
                    mars=mars,
                    outer_route=:none
                )
                seed_rows[i] = row
                if row.solve_success
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): total=$(round(row.total_time_s; digits=3)) s"
                else
                    seed_msgs[i] = "    seed $(i)/$(length(seeds))=$(seed): $(err)"
                end
            end
        end

        for i in eachindex(seeds)
            push!(rows, seed_rows[i])
            println(seed_msgs[i])
        end
        _record_outer_route_feedback!(warmup_case, seed_rows; route=mc_backend)
    end
    return nothing
end

function run_benchmarks(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    selected = spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
    selected = _split_rollout_benchmark_cases(selected)
    selected = _multirate_rollout_benchmark_cases(selected)
    rows = NamedTuple[]
    total = length(selected)
    backend = perf_parallel_backend()
    if backend == :auto
        case_rows = Vector{Vector{NamedTuple}}(undef, total)
        chosen_routes = fill(:none, total)
        process_tasks = Tuple{Int, BenchmarkCase}[]
        thread_indices = Int[]
        serial_indices = Int[]

        for (idx, case) in enumerate(selected)
            route = auto_backend_for_case(case; spec=spec)
            chosen_routes[idx] = route
            if route == :process
                push!(process_tasks, (idx, case))
            elseif route == :threads
                if _case_outer_threads_safe(case)
                    push!(thread_indices, idx)
                else
                    push!(serial_indices, idx)
                    chosen_routes[idx] = :none
                end
            else
                push!(serial_indices, idx)
            end
        end

        if !isempty(process_tasks)
            ensure_perf_workers!()
            process_rows = pmap(process_tasks) do task
                idx, case = task
                perf_worker_run_case_batch(case, spec, idx, total, :process)
            end
            for (k, task) in enumerate(process_tasks)
                idx = task[1]
                case_rows[idx] = process_rows[k]
            end
        end

        if !isempty(thread_indices)
            for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
                withenv(env_pairs...) do
                    Threads.@threads for j in eachindex(payload)
                        idx, plan = payload[j]
                        local_rows = NamedTuple[]
                        run_case_batch!(
                            local_rows,
                            selected[idx],
                            spec,
                            idx,
                            total;
                            outer_route=:threads,
                            plan=plan,
                            apply_env=false
                        )
                        case_rows[idx] = local_rows
                    end
                end
            end
        end

        for idx in serial_indices
            local_rows = NamedTuple[]
            run_case_batch!(local_rows, selected[idx], spec, idx, total; outer_route=:none)
            case_rows[idx] = local_rows
        end

        for idx in eachindex(case_rows)
            _record_outer_route_feedback!(selected[idx], case_rows[idx]; route=chosen_routes[idx])
            append!(rows, case_rows[idx])
        end
    elseif backend == :process
        ensure_perf_workers!()
        tasks = collect(enumerate(selected))
        case_rows = pmap(tasks) do task
            idx, case = task
            perf_worker_run_case_batch(case, spec, idx, total, :process)
        end
        for idx in eachindex(case_rows)
            _record_outer_route_feedback!(selected[idx], case_rows[idx]; route=:process)
            append!(rows, case_rows[idx])
        end
    elseif backend == :threads
        case_rows = Vector{Vector{NamedTuple}}(undef, total)
        thread_indices, serial_indices = _split_threadsafe_indices(selected, collect(eachindex(selected)))
        for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
            withenv(env_pairs...) do
                Threads.@threads for j in eachindex(payload)
                    idx, plan = payload[j]
                    local_rows = NamedTuple[]
                    run_case_batch!(
                        local_rows,
                        selected[idx],
                        spec,
                        idx,
                        total;
                        outer_route=:threads,
                        plan=plan,
                        apply_env=false
                    )
                    case_rows[idx] = local_rows
                end
            end
        end
        for idx in serial_indices
            local_rows = NamedTuple[]
            run_case_batch!(local_rows, selected[idx], spec, idx, total; outer_route=:none)
            case_rows[idx] = local_rows
        end
        for idx in eachindex(case_rows)
            route = idx in serial_indices ? :none : :threads
            _record_outer_route_feedback!(selected[idx], case_rows[idx]; route=route)
            append!(rows, case_rows[idx])
        end
    else
        for (idx, case) in enumerate(selected)
            local_rows = NamedTuple[]
            run_case_batch!(local_rows, case, spec, idx, total; outer_route=:none)
            _record_outer_route_feedback!(case, local_rows; route=:none)
            append!(rows, local_rows)
        end
    end

    run_montecarlo_batch!(rows, spec, planet)
    return DataFrame(rows)
end
@inline function selected_cases(spec::ProfileSpec, cases::Vector{BenchmarkCase})::Vector{BenchmarkCase}
    return spec.name == "full" ? cases : [c for c in cases if c.run_in_quick]
end

function measure_per_orbit_scenario(
    base_case::BenchmarkCase,
    spec::ProfileSpec,
    period_s::Float64,
    orbit_counts::Vector{Int};
    outer_route::Symbol=:none,
    apply_env::Bool=true
)
    rows = NamedTuple[]
    logs = String[]
    stream_logs = _perf_stream_orbit_logs()
    for orbit_count in orbit_counts
        mission_time = orbit_count * period_s
        args_template = deepcopy(base_case.args_template)
        args_template = SimulationConfiguration(
            file_paths=args_template.file_paths,
            simulation_settings=args_template.simulation_settings,
            mission_configuration=MissionConfiguration(
                mission_type=args_template.mission_configuration.mission_type,
                keplerian=args_template.mission_configuration.keplerian,
                number_of_orbits=args_template.mission_configuration.number_of_orbits,
                mission_time=mission_time,
                orientation_sim=args_template.mission_configuration.orientation_sim,
                num_steps_to_save=args_template.mission_configuration.num_steps_to_save
            ),
            environment_model=args_template.environment_model,
            dynamics_model=args_template.dynamics_model,
            guidance_model=args_template.guidance_model,
            navigation_model=args_template.navigation_model,
            control_model=args_template.control_model,
            initial_time=args_template.initial_time,
            integration_tolerances=args_template.integration_tolerances
        )

        case = BenchmarkCase(
            name=base_case.name,
            category=base_case.category,
            description=base_case.description,
            args_template=args_template,
            run_in_quick=base_case.run_in_quick,
            solver_mode_override=base_case.solver_mode_override,
            split_imex_solver_override=base_case.split_imex_solver_override,
            entry_target_count_override=base_case.entry_target_count_override
        )

        plan = parallel_priority_plan(case, outer_route)
        run_case = () -> begin
            if stream_logs
                println("    mission_time_multiplier=x$(orbit_count): warmup x$(spec.warmup), repeats x$(spec.repeats)")
                flush(stdout)
            end
            run_warmup(case, spec.warmup, spec.name)
            for rep in 1:spec.repeats
                last_row = nothing
                for attempt in 1:spec.max_attempts
                    row = measure_case(case, spec.name, rep; attempt=attempt, plan=plan)
                    row_orbit = merge(
                        row,
                        (
                            orbit_count=orbit_count,
                            mission_time_multiplier=orbit_count,
                            orbital_period_s=period_s
                        )
                    )
                    last_row = row_orbit
                    if row_orbit.solve_success
                        push!(rows, row_orbit)
                        line = "    mission_time_multiplier=x$(orbit_count) repeat $(rep)/$(spec.repeats): total=$(round(row_orbit.total_time_s; digits=3)) s"
                        push!(logs, line)
                        if stream_logs
                            println(line)
                            flush(stdout)
                        end
                        break
                    end
                end
                if !(last_row === nothing) && !last_row.solve_success
                    push!(rows, last_row)
                    line = "    mission_time_multiplier=x$(orbit_count) repeat $(rep)/$(spec.repeats): failed after $(spec.max_attempts) attempts, retcode=$(last_row.solve_retcode)"
                    push!(logs, line)
                    if stream_logs
                        println(line)
                        flush(stdout)
                    end
                end
            end
        end
        if apply_env
            env_pairs = parallel_priority_env_pairs(plan)
            withenv(env_pairs...) do
                run_case()
            end
        else
            run_case()
        end
    end
    return rows, logs
end

function run_montecarlo_per_orbit!(
    rows::Vector{NamedTuple},
    spec::ProfileSpec,
    planet::Earth,
    period_s::Float64,
    orbit_counts::Vector{Int}
)
    seeds = collect(1:spec.montecarlo_samples)
    scenarios = _active_montecarlo_scenarios()
    scenario_names = join([s.name for s in scenarios], ", ")
    println("  montecarlo scenarios (mission-time sweep, seeds=$(length(seeds))): $(scenario_names)")
    backend = perf_parallel_backend()
    mars = perf_worker_mars()
    workers_ready = false

    for scenario in scenarios
        variant = scenario.variant
        println("  scenario $(scenario.name)")
        for orbit_count in orbit_counts
            mission_time = orbit_count * period_s
            orbit_case = make_montecarlo_case(first(seeds), mission_time, variant, planet; mars=mars)
            mc_backend = backend == :auto ? auto_backend_for_case(orbit_case; spec=spec) : backend
            if mc_backend == :process && !workers_ready
                ensure_perf_workers!()
                workers_ready = true
            end
            println("    mission_time_multiplier=x$(orbit_count)")
            orbit_rows = Vector{NamedTuple}(undef, length(seeds))
            orbit_msgs = Vector{String}(undef, length(seeds))

            if mc_backend == :process
                seed_results = pmap(seed -> perf_worker_measure_montecarlo_seed(spec, mission_time, seed, variant, :process), seeds)
                for i in eachindex(seeds)
                    seed = seeds[i]
                    row, err = seed_results[i]
                    row_orbit = merge(row, (orbit_count=orbit_count, mission_time_multiplier=orbit_count, orbital_period_s=period_s))
                    orbit_rows[i] = row_orbit
                    if row_orbit.solve_success
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                    else
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                    end
                end
            elseif mc_backend == :threads
                threaded_plan = parallel_priority_plan(orbit_case, :threads)
                threaded_env = parallel_priority_env_pairs(threaded_plan)
                withenv(threaded_env...) do
                    Threads.@threads for i in eachindex(seeds)
                        seed = seeds[i]
                        row, err = measure_montecarlo_seed(
                            spec,
                            planet,
                            mission_time,
                            seed;
                            variant=variant,
                            mars=mars,
                            outer_route=:threads,
                            plan=threaded_plan,
                            apply_env=false
                        )
                        row_orbit = merge(row, (orbit_count=orbit_count, mission_time_multiplier=orbit_count, orbital_period_s=period_s))
                        orbit_rows[i] = row_orbit
                        if row_orbit.solve_success
                            orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                        else
                            orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                        end
                    end
                end
            else
                for i in eachindex(seeds)
                    seed = seeds[i]
                    row, err = measure_montecarlo_seed(
                        spec,
                        planet,
                        mission_time,
                        seed;
                        variant=variant,
                        mars=mars,
                        outer_route=:none
                    )
                    row_orbit = merge(row, (orbit_count=orbit_count, mission_time_multiplier=orbit_count, orbital_period_s=period_s))
                    orbit_rows[i] = row_orbit
                    if row_orbit.solve_success
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): total=$(round(row_orbit.total_time_s; digits=3)) s"
                    else
                        orbit_msgs[i] = "      seed $(i)/$(length(seeds))=$(seed): $(err)"
                    end
                end
            end

            for i in eachindex(seeds)
                push!(rows, orbit_rows[i])
                println(orbit_msgs[i])
            end
            _record_outer_route_feedback!(orbit_case, orbit_rows; route=mc_backend)
        end
    end
    return nothing
end
function run_per_orbit_for_scenarios(spec::ProfileSpec, cases::Vector{BenchmarkCase}, planet::Earth)::DataFrame
    baseline_sc = make_spacecraft(planet; id=1, with_panel=false)
    period_s = orbital_period_seconds(baseline_sc, planet)
    orbit_counts = spec.name == "full" ? collect(1:5) : collect(1:3)
    include_control_stress = _include_control_stress_per_orbit()
    selected_base = include_control_stress ? selected_cases(spec, cases) : [c for c in selected_cases(spec, cases) if c.category != "control_stress"]
    # Entry scenarios use entry-interface counting semantics, not baseline-period multipliers.
    selected = [c for c in selected_base if c.category != "entry"]
    excluded_entry = length(selected_base) - length(selected)

    println(
        "[mission-time-sweep] scenarios=$(length(selected)), baseline period=$(round(period_s; digits=3)) s, " *
        "multipliers=x$(first(orbit_counts)):x$(last(orbit_counts)), include_control_stress=$(include_control_stress), " *
        "entry_excluded=$(excluded_entry)"
    )
    rows = NamedTuple[]
    backend = perf_parallel_backend()

    if backend == :auto
        scenario_rows = Vector{Vector{NamedTuple}}(undef, length(selected))
        scenario_logs = Vector{Vector{String}}(undef, length(selected))
        chosen_routes = fill(:none, length(selected))
        process_tasks = Tuple{Int, BenchmarkCase}[]
        thread_indices = Int[]
        serial_indices = Int[]

        for (idx, base_case) in enumerate(selected)
            route = auto_backend_for_case(base_case; spec=spec)
            chosen_routes[idx] = route
            if route == :process
                push!(process_tasks, (idx, base_case))
            elseif route == :threads
                if _case_outer_threads_safe(base_case)
                    push!(thread_indices, idx)
                else
                    push!(serial_indices, idx)
                    chosen_routes[idx] = :none
                end
            else
                push!(serial_indices, idx)
            end
        end

        if !isempty(process_tasks)
            ensure_perf_workers!()
            process_results = pmap(process_tasks) do task
                idx, base_case = task
                local_rows, local_logs = perf_worker_measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts, :process)
                return (rows=local_rows, logs=local_logs)
            end
            for (k, task) in enumerate(process_tasks)
                idx = task[1]
                scenario_rows[idx] = process_results[k].rows
                scenario_logs[idx] = process_results[k].logs
            end
        end

        if !isempty(thread_indices)
            for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
                withenv(env_pairs...) do
                    Threads.@threads for j in eachindex(payload)
                        idx, _ = payload[j]
                        local_rows, local_logs = measure_per_orbit_scenario(
                            selected[idx],
                            spec,
                            period_s,
                            orbit_counts;
                            outer_route=:threads,
                            apply_env=false
                        )
                        scenario_rows[idx] = local_rows
                        scenario_logs[idx] = local_logs
                    end
                end
            end
        end

        for idx in serial_indices
            local_rows, local_logs = measure_per_orbit_scenario(selected[idx], spec, period_s, orbit_counts; outer_route=:none)
            scenario_rows[idx] = local_rows
            scenario_logs[idx] = local_logs
        end

        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            _record_outer_route_feedback!(base_case, scenario_rows[idx]; route=chosen_routes[idx])
            append!(rows, scenario_rows[idx])
            for log_line in scenario_logs[idx]
                println(log_line)
            end
        end
    elseif backend == :process
        ensure_perf_workers!()
        tasks = collect(enumerate(selected))
        scenario_results = pmap(tasks) do task
            idx, base_case = task
            local_rows, local_logs = perf_worker_measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts, :process)
            return (rows=local_rows, logs=local_logs)
        end
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            _record_outer_route_feedback!(base_case, scenario_results[idx].rows; route=:process)
            append!(rows, scenario_results[idx].rows)
            for log_line in scenario_results[idx].logs
                println(log_line)
            end
        end
    elseif backend == :threads
        scenario_rows = Vector{Vector{NamedTuple}}(undef, length(selected))
        scenario_logs = Vector{Vector{String}}(undef, length(selected))
        thread_indices, serial_indices = _split_threadsafe_indices(selected, collect(eachindex(selected)))
        for (env_pairs, payload) in _thread_plan_groups(selected, thread_indices, :threads)
            withenv(env_pairs...) do
                Threads.@threads for j in eachindex(payload)
                    idx, _ = payload[j]
                    local_rows, local_logs = measure_per_orbit_scenario(
                        selected[idx],
                        spec,
                        period_s,
                        orbit_counts;
                        outer_route=:threads,
                        apply_env=false
                    )
                    scenario_rows[idx] = local_rows
                    scenario_logs[idx] = local_logs
                end
            end
        end
        for idx in serial_indices
            local_rows, local_logs = measure_per_orbit_scenario(
                selected[idx],
                spec,
                period_s,
                orbit_counts;
                outer_route=:none
            )
            scenario_rows[idx] = local_rows
            scenario_logs[idx] = local_logs
        end
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            route = idx in serial_indices ? :none : :threads
            _record_outer_route_feedback!(base_case, scenario_rows[idx]; route=route)
            append!(rows, scenario_rows[idx])
            for log_line in scenario_logs[idx]
                println(log_line)
            end
        end
    else
        for (idx, base_case) in enumerate(selected)
            println("  scenario $(idx)/$(length(selected)) = $(base_case.name)")
            local_rows, local_logs = measure_per_orbit_scenario(base_case, spec, period_s, orbit_counts; outer_route=:none)
            _record_outer_route_feedback!(base_case, local_rows; route=:none)
            append!(rows, local_rows)
            for log_line in local_logs
                println(log_line)
            end
        end
    end

    run_montecarlo_per_orbit!(rows, spec, planet, period_s, orbit_counts)
    return DataFrame(rows)
end
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

@inline function _summary_group_keys(df::DataFrame, base_keys::Vector{Symbol})::Vector{Symbol}
    keys = copy(base_keys)
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
    df_names = Set(Symbol.(names(df)))
    for key in optional_keys
        if key in df_names
            push!(keys, key)
        end
    end
    return keys
end

function summarize_per_orbit_results(orbit_raw_df::DataFrame)::DataFrame
    sweep_multiplier_key = :mission_time_multiplier in names(orbit_raw_df) ? :mission_time_multiplier : :orbit_count
    keys = _summary_group_keys(
        orbit_raw_df,
        [:category, :scenario, :description, sweep_multiplier_key, :orbital_period_s, :dt_max_orbit_s]
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
        summary = leftjoin(counts, success_summary, on=keys)
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
            :density_parallel_mode,
            :control_parallel_mode,
            :multibody_parallel_mode,
            :dt_max_orbit_s,
            :dynamic_effectors,
            :control_effectors
        ]
    )
    counts = combine(
        groupby(raw_df, keys),
        nrow => :samples_total,
        :solve_success => (v -> count(identity, v)) => :samples_success
    )
    counts[!, :samples_failed] = counts.samples_total .- counts.samples_success
    counts[!, :success_rate] = Float64.(counts.samples_success) ./ Float64.(counts.samples_total)

    success_df = raw_df[raw_df.solve_success .== true, :]
    metric_cols = [
        :samples,
        :copy_time_mean_s,
        :solve_time_mean_s,
        :total_time_mean_s,
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

    for t in sample_times
        u_ref = _solution_state_at(sol_ref, t)
        u_cmp = _solution_state_at(sol_cmp, t)
        for sat_idx in 1:n_sats
            sc_ref = u_ref.sc[sat_idx]
            sc_cmp = u_cmp.sc[sat_idx]
            pos_rel_max = max(pos_rel_max, _relative_vector_delta(sc_ref.pos, sc_cmp.pos))
            vel_rel_max = max(vel_rel_max, _relative_vector_delta(sc_ref.vel, sc_cmp.vel))
            if orientation
                q_angle_max_rad = max(q_angle_max_rad, _quaternion_angle_delta_rad(sc_ref.q, sc_cmp.q))
                omega_rel_max = max(omega_rel_max, _relative_vector_delta(sc_ref.ω, sc_cmp.ω))
            end
        end
    end

    return (
        t_start=t_start,
        t_end=t_end,
        sample_count=sample_count,
        pos_rel_max=pos_rel_max,
        vel_rel_max=vel_rel_max,
        q_angle_max_rad=orientation ? q_angle_max_rad : missing,
        omega_rel_max=orientation ? omega_rel_max : missing
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
    orbit_summary_df::DataFrame
)::Vector{String}
    !_plot_ready() && return String[]
    _ensure_runtime_plot_theme!()

    plot_artifacts = String[]
    success_summary = summary_df[summary_df.samples_success .> 0, :]

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

    # 3) Runtime variability by scenario.
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

    # 4) Configuration copy + solve breakdown.
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

    # 5) Allocation footprint by scenario (memory + call count).
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

    # 6) Integrator workload and rejection pressure.
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

    # 7) Throughput ranking with per-satellite markers.
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

    # 8) Satellite-count scaling (measured vs ideal linear).
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

    # 9) Dynamics fidelity ladder (absolute + relative).
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

    # 10) Monte Carlo runtime distribution with mean and p90 lines.
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

    # 11) Monte Carlo seed trace with mean and p90 lines.
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

    # 12) Mission-time sweep runtime scaling.
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

    # 13) Mission-time sweep efficiency scaling.
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

    # 14) Mission-time sweep time heatmap.
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

    # 15) SPICE call budget by scenario (N-body + SRP + planet frame).
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
    orbit_summary_df::DataFrame;
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
    density_backend_breakdown_csv_path::Union{Nothing, String}=nothing
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
                "mission_time_sweep=`$(_fmt(orbit_stage_s))`, total=`$(_fmt(total_stage_s))`"
            )
        end
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
    end
end

function main()
    spec, outdir = parse_cli()
    mkpath(outdir)
    outer_route_state = _load_outer_route_history!(spec, outdir)

    solver_mode_default = _perf_default_solver_mode(spec.name)
    solver_mode_effective = _perf_solver_mode_env(spec.name)
    sat_threshold = _priority_inner_sat_threshold()
    link_threshold = _priority_inner_link_threshold()
    outer_light_sat = _priority_outer_light_sat_threshold()
    outer_light_link = _priority_outer_light_link_threshold()
    outer_light_mission_s = _priority_outer_light_mission_threshold_s()

    println("Performance runtime analysis profile: $(spec.name)")
    println("Output directory: $outdir")
    println("Solver mode default: $(solver_mode_default)")
    println("Solver mode effective: $(solver_mode_effective)")
    println(
        "Outer-route adaptive=$(_outer_route_adaptive_enabled() ? "on" : "off"), " *
        "min_samples=$(_outer_route_min_samples()), " *
        "mc_process_min_samples=$(_outer_route_mc_process_min_samples()), " *
        "mc_process_min_mission_s=$(round(_outer_route_mc_process_min_mission_s(); digits=3))"
    )
    if outer_route_state.persist
        if outer_route_state.reset_requested
            println("Outer-route state cache reset requested; starting from empty history.")
        elseif outer_route_state.loaded_rows > 0
            println(
                "Outer-route state cache loaded rows=$(outer_route_state.loaded_rows), " *
                "signatures=$(outer_route_state.loaded_signatures) from $(outer_route_state.path)"
            )
        else
            println("Outer-route state cache path: $(outer_route_state.path) (no prior state loaded)")
        end
    else
        println("Outer-route state cache persistence: off")
    end
    println(
        "Priority thresholds: inner_sat=$(sat_threshold), inner_link=$(link_threshold), " *
        "outer_light_sat=$(outer_light_sat), outer_light_link=$(outer_light_link), " *
        "outer_light_mission_s=$(round(outer_light_mission_s; digits=3))"
    )

    planet = Earth("", SPICE_PATH)
    cases = build_cases(spec, planet)
    bench_started_ns = time_ns()
    raw_df = run_benchmarks(spec, cases, planet)
    summary_df = summarize_results(raw_df)
    density_backend_breakdown_df = summarize_density_backend_breakdown(raw_df)
    bench_elapsed_s = (time_ns() - bench_started_ns) / 1e9

    split_gate_df = nothing
    split_gate_csv_path = nothing
    split_gate_report_path = nothing
    split_gate_elapsed_s = 0.0
    if _split_rollout_enabled()
        split_gate_started_ns = time_ns()
        split_gate_result = evaluate_split_rollout_gate(spec, cases, outdir)
        split_gate_elapsed_s = (time_ns() - split_gate_started_ns) / 1e9
        split_gate_df = split_gate_result.df
        split_gate_csv_path = split_gate_result.csv_path
        split_gate_report_path = split_gate_result.report_path
    end

    multirate_gate_df = nothing
    multirate_gate_csv_path = nothing
    multirate_gate_report_path = nothing
    multirate_gate_elapsed_s = 0.0
    if _multirate_rollout_enabled()
        multirate_gate_started_ns = time_ns()
        multirate_gate_result = evaluate_multirate_rollout_gate(spec, cases, outdir)
        multirate_gate_elapsed_s = (time_ns() - multirate_gate_started_ns) / 1e9
        multirate_gate_df = multirate_gate_result.df
        multirate_gate_csv_path = multirate_gate_result.csv_path
        multirate_gate_report_path = multirate_gate_result.report_path
    end

    orbit_started_ns = time_ns()
    orbit_raw_df = run_per_orbit_for_scenarios(spec, cases, planet)
    orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)
    orbit_elapsed_s = (time_ns() - orbit_started_ns) / 1e9
    total_elapsed_s = bench_elapsed_s + split_gate_elapsed_s + multirate_gate_elapsed_s + orbit_elapsed_s
    outer_route_state_saved = _save_outer_route_history!(spec, outer_route_state.path, outer_route_state.persist)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    raw_path = joinpath(outdir, "runtime_raw_$(spec.name)_$(stamp).csv")
    summary_path = joinpath(outdir, "runtime_summary_$(spec.name)_$(stamp).csv")
    orbit_raw_path = joinpath(outdir, "runtime_per_orbit_raw_$(spec.name)_$(stamp).csv")
    orbit_summary_path = joinpath(outdir, "runtime_per_orbit_summary_$(spec.name)_$(stamp).csv")
    stage_timing_path = joinpath(outdir, "runtime_stage_timing_$(spec.name)_$(stamp).csv")
    hardware_info_path = joinpath(outdir, "runtime_hardware_info_$(spec.name)_$(stamp).csv")
    inner_hint_layer_path = joinpath(outdir, "runtime_inner_hint_layers_$(spec.name)_$(stamp).csv")
    density_backend_breakdown_path = joinpath(outdir, "runtime_density_backend_breakdown_$(spec.name)_$(stamp).csv")
    report_path = joinpath(outdir, "runtime_report_$(spec.name)_$(stamp).md")
    plot_paths = generate_runtime_plots(outdir, spec, stamp, raw_df, summary_df, orbit_summary_df)

    stage_names = ["run_benchmarks"]
    stage_elapsed = [bench_elapsed_s]
    if _split_rollout_enabled()
        push!(stage_names, "run_split_rollout_gate")
        push!(stage_elapsed, split_gate_elapsed_s)
    end
    if _multirate_rollout_enabled()
        push!(stage_names, "run_multirate_rollout_gate")
        push!(stage_elapsed, multirate_gate_elapsed_s)
    end
    push!(stage_names, "run_per_orbit")
    push!(stage_names, "total")
    push!(stage_elapsed, orbit_elapsed_s)
    push!(stage_elapsed, total_elapsed_s)
    stage_timing_df = DataFrame(stage=stage_names, elapsed_s=stage_elapsed)
    hw = _runtime_hardware_snapshot()
    hardware_info_df = DataFrame([
        (
            profile=spec.name,
            machine_label=hw.machine_label,
            hardware_class=hw.hardware_class,
            host_name=hw.host_name,
            cpu_name=hw.cpu_name,
            cpu_threads=hw.cpu_threads,
            julia_threads=hw.julia_threads,
            os=hw.os,
            arch=hw.arch
        )
    ])
    inner_hint_layer_df = _inner_hint_layer_report_df(spec, hw)

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    CSV.write(orbit_raw_path, orbit_raw_df)
    CSV.write(orbit_summary_path, orbit_summary_df)
    CSV.write(stage_timing_path, stage_timing_df)
    CSV.write(hardware_info_path, hardware_info_df)
    CSV.write(inner_hint_layer_path, inner_hint_layer_df)
    CSV.write(density_backend_breakdown_path, density_backend_breakdown_df)
    write_report(
        report_path,
        spec,
        raw_df,
        summary_df,
        orbit_summary_df;
        plot_paths=plot_paths,
        split_gate_df=split_gate_df,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_report_path=split_gate_report_path,
        multirate_gate_df=multirate_gate_df,
        multirate_gate_csv_path=multirate_gate_csv_path,
        multirate_gate_report_path=multirate_gate_report_path,
        stage_timing_df=stage_timing_df,
        inner_hint_layer_df=inner_hint_layer_df,
        inner_hint_layer_csv_path=inner_hint_layer_path,
        density_backend_breakdown_df=density_backend_breakdown_df,
        density_backend_breakdown_csv_path=density_backend_breakdown_path
    )

    println("Analysis complete.")
    println("Raw results: $raw_path")
    println("Summary: $summary_path")
    println("Mission-time-sweep raw: $orbit_raw_path")
    println("Mission-time-sweep summary: $orbit_summary_path")
    println("Stage timing: $stage_timing_path")
    println("Hardware info: $hardware_info_path")
    println("Inner hint layers: $inner_hint_layer_path")
    println("Density backend breakdown: $density_backend_breakdown_path")
    if !(split_gate_df === nothing)
        pass_count = (:pass_all in names(split_gate_df)) ? count(Bool.(split_gate_df.pass_all)) : 0
        println("Split rollout gate: $(pass_count)/$(nrow(split_gate_df)) pass")
    end
    if !(split_gate_csv_path === nothing)
        println("Split rollout gate CSV: $(split_gate_csv_path)")
    end
    if !(split_gate_report_path === nothing)
        println("Split rollout gate report: $(split_gate_report_path)")
    end
    if !(multirate_gate_df === nothing)
        pass_count = (:pass_all in names(multirate_gate_df)) ? count(Bool.(multirate_gate_df.pass_all)) : 0
        println("Multirate rollout gate: $(pass_count)/$(nrow(multirate_gate_df)) pass")
    end
    if !(multirate_gate_csv_path === nothing)
        println("Multirate rollout gate CSV: $(multirate_gate_csv_path)")
    end
    if !(multirate_gate_report_path === nothing)
        println("Multirate rollout gate report: $(multirate_gate_report_path)")
    end
    if outer_route_state.persist
        println(
            "Outer-route state cache: $(outer_route_state_saved.path) " *
            "(rows_saved=$(outer_route_state_saved.rows), signatures=$(outer_route_state_saved.signatures))"
        )
    end
    println("Plots generated: $(length(plot_paths))")
    println("Report: $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
