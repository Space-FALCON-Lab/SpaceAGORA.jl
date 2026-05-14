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

@inline function _mode_token(v)::String
    return lowercase(strip(string(v)))
end

@inline function _is_adaptive_mode_token(v)::Bool
    token = _mode_token(v)
    return token in ("outer_inner_adaptive", "outer_inner_full_smart")
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

@inline function _perf_smoke_mode()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PERF_SMOKE", "0")))
    return raw in ("1", "true", "yes", "on")
end

function _perf_case_name_filter()::Vector{String}
    raw = strip(get(ENV, "SPACEAGORA_PERF_CASES", ""))
    isempty(raw) && return String[]
    names = String[]
    seen = Set{String}()
    for part in split(raw, ",")
        name = strip(part)
        isempty(name) && continue
        if !(name in seen)
            push!(names, name)
            push!(seen, name)
        end
    end
    return names
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
    entry_target_count_override::Union{Nothing, Int}=nothing,
    solver_cache=nothing
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
            return_solver_metadata=return_solver_metadata,
            solver_cache=solver_cache
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
        spec = ProfileSpec(
            name="full",
            repeats=5,
            warmup=2,
            max_attempts=4,
            mission_short_s=3600.0,
            mission_long_s=14400.0,
            montecarlo_samples=50,
            montecarlo_mission_s=3600.0
        )
        if _perf_smoke_mode()
            return ProfileSpec(
                name=spec.name,
                repeats=1,
                warmup=0,
                max_attempts=2,
                mission_short_s=spec.mission_short_s,
                mission_long_s=spec.mission_long_s,
                montecarlo_samples=2
                ,
                montecarlo_mission_s=spec.montecarlo_mission_s
            )
        end
        return spec
    elseif name == "quick"
        spec = ProfileSpec(
            name="quick",
            repeats=3,
            warmup=1,
            max_attempts=4,
            mission_short_s=1800.0,
            mission_long_s=7200.0,
            montecarlo_samples=50,
            montecarlo_mission_s=1800.0
        )
        if _perf_smoke_mode()
            return ProfileSpec(
                name=spec.name,
                repeats=1,
                warmup=0,
                max_attempts=2,
                mission_short_s=spec.mission_short_s,
                mission_long_s=spec.mission_long_s,
                montecarlo_samples=2
                ,
                montecarlo_mission_s=spec.montecarlo_mission_s
            )
        end
        return spec
    else
        throw(ArgumentError("Unsupported profile '$name'. Valid values: quick, full."))
    end
end

@inline function perf_parallel_enabled()::Bool
    mode = lowercase(strip(get(ENV, "SPACEAGORA_PERF_PARALLEL", "auto")))
    return !(mode in ("0", "false", "off", "no"))
end

@inline perf_threads_available()::Bool = Threads.nthreads() > 1


@inline function perf_parallel_backend()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_PERF_PARALLEL_BACKEND", "auto")))
    if mode == "auto" && !isempty(_perf_case_name_filter())
        return :none
    end
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

@inline function _exclude_entry_scenarios()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_EXCLUDE_ENTRY_SCENARIOS", false)
end

@inline function _include_montecarlo_scenarios()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_INCLUDE_MONTECARLO", true)
end

@inline function _include_mission_time_sweep()::Bool
    return _parse_bool_env("SPACEAGORA_PERF_INCLUDE_MISSION_TIME_SWEEP", true)
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
