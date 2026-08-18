@inline function parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function parse_parallel_mode_env(name::String; default::String="auto")::Symbol
    mode = lowercase(strip(get(ENV, name, default)))
    if mode in ("off", "none", "serial", "0", "false", "no")
        return :off
    elseif mode in ("on", "thread", "threads", "1", "true", "yes")
        return :on
    elseif mode == "auto"
        return :auto
    end
    throw(ArgumentError("Unsupported $name='$mode'. Use one of: off, auto, on."))
end

@inline function parse_thread_threshold_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    threshold = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    return max(1, threshold)
end

@inline function parse_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a float, got '$raw'"))
    end
    return value
end

@inline function parse_nonnegative_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    return max(0, value)
end

@inline function inner_scheduler_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER", "static")))
    if mode in ("static", "strided")
        return :static
    elseif mode == "dynamic"
        return :dynamic
    end
    throw(ArgumentError("Unsupported SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER='$mode'. Use one of: static, dynamic."))
end

@inline function inner_dynamic_chunk_size()::Int
    return parse_thread_threshold_env("SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK", 1)
end

@inline function _safe_token(raw)::String
    token = lowercase(strip(String(raw)))
    token = replace(token, r"[^a-z0-9._-]+" => "_")
    return isempty(token) ? "default" : token
end

@inline persistent_hints_enabled()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS", false)

@inline function persistent_hints_persist_enabled()::Bool
    return parse_bool_env("SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST", persistent_hints_enabled())
end

@inline function persistent_hints_state_reset_requested()::Bool
    return parse_bool_env("SPACEAGORA_PARALLEL_POLICY_STATE_RESET", false)
end

@inline function persistent_hints_exploration()::Float64
    c = parse_float_env("SPACEAGORA_PARALLEL_POLICY_HINT_EXPLORATION", 1.5)
    return c > 0.0 ? c : 1.5
end

@inline function persistent_hints_min_samples()::Int
    return parse_thread_threshold_env("SPACEAGORA_PARALLEL_POLICY_HINT_MIN_SAMPLES", 2)
end

@inline function _persistent_hint_default_path()::String
    profile = _safe_token(get(ENV, "SPACEAGORA_PARALLEL_PROFILE", "default"))
    machine = _safe_token(get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", "default"))
    threads = Threads.nthreads()
    return joinpath(
        pwd(),
        "output",
        "parallel_policy_state",
        "inner_policy_state_$(profile)_$(machine)_t$(threads).toml"
    )
end

@inline function _persistent_hint_path()::String
    override = strip(get(ENV, "SPACEAGORA_PARALLEL_POLICY_STATE_PATH", ""))
    if isempty(override)
        return normpath(_persistent_hint_default_path())
    end
    return normpath(isabspath(override) ? override : joinpath(pwd(), override))
end

@inline outer_parallel_active()::Bool = parse_bool_env("SPACEAGORA_OUTER_PARALLEL_ACTIVE", false)

@inline function _default_thread_pool_size()::Int
    try
        return Threads.nthreads()
    catch
        return Threads.nthreads()
    end
end

@inline function effective_inner_thread_budget()::Int
    raw = strip(get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", "0"))
    budget = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_INNER_THREAD_BUDGET must be an integer, got '$raw'"))
    end
    available = _default_thread_pool_size()
    if budget <= 0
        return max(1, available)
    end
    return max(1, min(available, budget))
end

@inline function auto_thread_min_budget()::Int
    return parse_thread_threshold_env("SPACEAGORA_AUTO_THREAD_MIN_BUDGET", 4)
end

@inline function auto_thread_min_budget(source::Symbol)::Int
    default_budget = auto_thread_min_budget()
    if source == :density_callback
        return parse_thread_threshold_env(
            "SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET",
            max(default_budget, 16),
        )
    elseif source == :density_callback_lockfree
        # A lock-free density model (e.g. GRAMAtmosphereModelSurrogate) has no
        # global-lock oversubscription cost to guard against -- the 16-thread
        # floor above exists for native/locked GRAM specifically, and gates
        # the lock-free surrogate for no reason. Falls through to the general
        # default (4) like :multibody/:dynamic_effectors/:control_callback.
        return parse_thread_threshold_env(
            "SPACEAGORA_DENSITY_CALLBACK_LOCKFREE_AUTO_THREAD_MIN_BUDGET",
            default_budget,
        )
    elseif source == :thermal_callback
        return parse_thread_threshold_env(
            "SPACEAGORA_THERMAL_CALLBACK_AUTO_THREAD_MIN_BUDGET",
            max(default_budget, 16),
        )
    end
    return default_budget
end

@inline function _telemetry_bucket(source::Symbol)::Symbol
    if source == :density_callback || source == :density_callback_lockfree
        return :density
    elseif source == :control_callback
        return :control
    elseif source == :multibody
        return :multibody
    end
    return :other
end

"""
    snapshot_policy_decision_env() -> PolicyDecisionEnvConfig

Resolve the env-derived knobs consulted on every `thread_policy_decision` call
into a typed snapshot.  Built once at run_simulation setup so hot paths avoid
process-global ENV access; values reflect ENV (and any active engine
overrides) at snapshot time.
"""
function snapshot_policy_decision_env()::PolicyDecisionEnvConfig
    return PolicyDecisionEnvConfig(
        effective_inner_thread_budget(),
        outer_parallel_active(),
        auto_thread_min_budget(),
        auto_thread_min_budget(:density_callback),
        auto_thread_min_budget(:density_callback_lockfree),
        auto_thread_min_budget(:thermal_callback),
        adaptive_policy_enabled(),
    )
end

@inline function _snapshot_auto_min_budget(env::PolicyDecisionEnvConfig, source::Symbol)::Int
    if source == :density_callback
        return env.auto_min_budget_density
    elseif source == :density_callback_lockfree
        return env.auto_min_budget_density_lockfree
    elseif source == :thermal_callback
        return env.auto_min_budget_thermal
    end
    return env.auto_min_budget_default
end

@inline adaptive_policy_enabled()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_ADAPTIVE", false)
@inline adaptive_window_size()::Int = parse_thread_threshold_env("SPACEAGORA_PARALLEL_POLICY_WINDOW", 8)
@inline adaptive_trim_quanta_budget()::Int = parse_nonnegative_int_env("SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA", 0)
@inline adaptive_bootstrap_threads()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS", true)
@inline adaptive_control_tail_guard()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD", false)
@inline adaptive_measured_reward_enabled()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD", persistent_hints_enabled())

@inline function adaptive_delta()::Float64
    δ = parse_float_env("SPACEAGORA_PARALLEL_POLICY_DELTA", 0.85)
    if !(0.0 < δ <= 1.0)
        throw(ArgumentError("SPACEAGORA_PARALLEL_POLICY_DELTA must satisfy 0 < δ <= 1, got '$δ'"))
    end
    return δ
end

@inline function adaptive_rho()::Float64
    ρ = parse_float_env("SPACEAGORA_PARALLEL_POLICY_RHO", 1.5)
    if !(ρ > 1.0)
        throw(ArgumentError("SPACEAGORA_PARALLEL_POLICY_RHO must satisfy ρ > 1, got '$ρ'"))
    end
    return ρ
end

@inline function _adaptive_desire_cap(pool_size::Int, ρ::Float64)::Int
    return max(1, ceil(Int, ρ * max(1, pool_size)))
end
