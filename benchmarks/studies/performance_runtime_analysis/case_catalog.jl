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
    elapsed_success_sq_sum_s = 0.0
    for row in rows
        if !hasproperty(row, :solve_success)
            continue
        end
        if row.solve_success === true
            if hasproperty(row, :total_time_s) && row.total_time_s isa Real && isfinite(Float64(row.total_time_s))
                elapsed_s = Float64(row.total_time_s)
                successes += 1
                elapsed_success_s += elapsed_s
                elapsed_success_sq_sum_s += elapsed_s^2
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
        elapsed_success_sq_sum_s=elapsed_success_sq_sum_s,
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
    serial_inner_off = plan.density_mode == "off" && plan.control_mode == "off" && plan.multibody_mode == "off"
    rhs_batch_mode = (plan.outer_route == :none && serial_inner_off) ? "off" : "auto"
    return Pair{String, String}[
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _existing_or_policy_env("SPACEAGORA_OUTER_PARALLEL_ACTIVE", outer_flag),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL", plan.density_mode),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => _existing_or_policy_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL", plan.control_mode),
        "SPACEAGORA_MULTIBODY_PARALLEL" => _existing_or_policy_env("SPACEAGORA_MULTIBODY_PARALLEL", plan.multibody_mode),
        "SPACEAGORA_RHS_BATCH_PARALLEL" => _existing_or_policy_env("SPACEAGORA_RHS_BATCH_PARALLEL", rhs_batch_mode),
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
    sc_thermal_aerobrake = [make_spacecraft(
        mars;
        id=1,
        with_panel=true,
        panel_count=16,
        ra_alt_m=4500e3,
        rp_alt_m=120e3,
        i_deg=93.0,
        ω_deg=80.0,
        Ω_deg=30.0,
        ν_deg=180.0,
        root_mass=340.0,
        root_area=5.0,
        panel_mass=6.0,
        panel_area=1.8,
        panel_offset_y=1.2
    )]
    sc_srp_heavy = [make_spacecraft(
        planet;
        id=1,
        with_panel=true,
        panel_count=8,
        ra_alt_m=700e3,
        rp_alt_m=680e3,
        i_deg=97.0,
        ω_deg=10.0,
        Ω_deg=45.0,
        ν_deg=30.0,
        root_mass=600.0,
        root_area=40.0,
        panel_mass=20.0,
        panel_area=35.0,
        panel_offset_y=2.5
    )]
    srp_heavy_effectors = (
        InverseSquaredGravityModel(),
        SolarRadiationPressureModel(1.2, 120.0),
        SolarRadiationPressureModel(1.8, 220.0)
    )
    sc_articulated_heavy = [make_spacecraft(
        planet;
        id=1,
        with_panel=true,
        panel_count=28,
        orientation_state=(q0, w0),
        ra_alt_m=520e3,
        rp_alt_m=500e3,
        ν_deg=160.0,
        root_mass=450.0,
        root_area=10.0,
        panel_mass=8.0,
        panel_area=3.0,
        panel_offset_y=2.0
    )]
    sc_multi_sat_control = make_constellation(planet, 8; with_panel=false)
    constellation_thruster = BaseThrusterModel(
        thrust=fill(0.18, 8),
        direction=fill(0.0, 8),
        Δv=fill(0.0, 8),
        start_burn_time=repeat([120.0, 540.0, 960.0, 1380.0], 2),
        stop_burn_time=repeat([180.0, 600.0, 1020.0, 1440.0], 2),
        Isp=fill(285.0, 8)
    )
    sc_long_constellation = make_constellation(planet, 12; with_panel=false)
    sc_effector_stress6 = make_constellation(planet, 6; with_panel=false)
    sc_effector_stress12 = make_constellation(planet, 12; with_panel=false)
    effector_stress_effectors = (
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
            name="thermal_aerobrake_mars_panel16",
            category="thermal_entry",
            description="1 articulated spacecraft (17 links) in Mars aerobraking regime (10 orbits, aero + MarsGRAM point-to-point) to stress thermal callback under entry-like heating",
            args_template=build_config(
                planet=mars,
                spacecraft=sc_thermal_aerobrake,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                mission_type=MissionOrbits,
                mission_keplerian=false,
                mission_orbits=10,
                dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
                density_model=deepcopy(mars_gram_point_density),
                dt_max_orbit=1.0
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="srp_heavy_high_area",
            category="srp_heavy",
            description="1 high-area spacecraft (9 links) with stacked SRP effectors to stress SRP-heavy workloads",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_srp_heavy,
                mission_time_s=spec.mission_short_s,
                orientation_sim=true,
                dynamic_effectors=srp_heavy_effectors,
                dt_max_orbit=1.0
            )
        ),
        BenchmarkCase(
            name="articulated_1sat_panel28_fullstack",
            category="articulated_multibody",
            description="1 heavily articulated spacecraft (29 links), orientation dynamics on, harmonics + aero active to stress multibody/link-level kernels",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_articulated_heavy,
                mission_time_s=min(spec.mission_short_s, 2400.0),
                orientation_sim=true,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20, AerodynamicCoefficientfM()),
                density_model=deepcopy(thermal_stress_density),
                dt_max_orbit=0.5
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="multi_sat_control_8sat_thruster",
            category="multi_sat_control",
            description="8-spacecraft constellation with active thruster control callbacks to capture multi-satellite + active control behavior",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_multi_sat_control,
                mission_time_s=min(spec.mission_short_s, 2400.0),
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredJ2GravityModel(),),
                control_effectors=(constellation_thruster,),
                control_rates=[0.2],
                dt_max_orbit=0.5
            )
        ),
        BenchmarkCase(
            name="long_constellation_12sat",
            category="long_constellation",
            description="12-spacecraft long-duration constellation with L20 harmonics + SRP and GRAM surrogate density",
            args_template=build_config(
                planet=planet,
                spacecraft=sc_long_constellation,
                mission_time_s=spec.mission_long_s,
                orientation_sim=false,
                dynamic_effectors=(InverseSquaredGravityModel(), harmonics20, SolarRadiationPressureModel(1.2, 12.0)),
                density_model=deepcopy(earth_gram_surrogate_density),
                dt_max_orbit=2.0
            ),
            run_in_quick=false
        ),
        BenchmarkCase(
            name="effector_6sat_dual_srp_stack",
            category="effector_stress",
            description="6 spacecraft with J2 + dual SRP effectors (no atmosphere/control) to stress dynamic effector reduction with outer routing off (r2)",
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
            description="12 spacecraft with J2 + dual SRP effectors (no atmosphere/control) for larger multi-satellite effector scaling",
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

