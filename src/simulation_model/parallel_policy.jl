module ParallelPolicy

using Base.Threads

export parse_bool_env, parse_parallel_mode_env, parse_thread_threshold_env
export outer_parallel_active, effective_inner_thread_budget, use_threads_policy
export thread_policy_decision, threaded_foreach, with_policy_context
export reset_policy_telemetry!, policy_telemetry_snapshot, record_policy_observation!

Base.@kwdef mutable struct AdaptiveControllerState
    desire::Int64 = 1
    window_calls::Int64 = 0
    window_allotment_sum::Int64 = 0
    window_useful_sum::Float64 = 0.0
    window_deprived_calls::Int64 = 0
    last_classification::Symbol = :none
    last_utilization::Float64 = 1.0
end

Base.@kwdef mutable struct PolicyTelemetry
    decisions_total::Int64 = 0
    threads_enabled_total::Int64 = 0
    adaptive_decisions_total::Int64 = 0
    density_decisions::Int64 = 0
    density_threads_enabled::Int64 = 0
    control_decisions::Int64 = 0
    control_threads_enabled::Int64 = 0
    multibody_decisions::Int64 = 0
    multibody_threads_enabled::Int64 = 0
    other_decisions::Int64 = 0
    other_threads_enabled::Int64 = 0
    last_source::Symbol = :none
    last_mode::Symbol = :none
    last_threshold::Int64 = 0
    last_num_items::Int64 = 0
    last_budget::Int64 = 0
    last_outer_active::Bool = false
    last_allow_with_outer::Bool = true
    last_heavy_only::Bool = false
    last_heavy_work::Bool = true
    last_use_threads::Bool = false
    last_adaptive_enabled::Bool = false
    last_desire::Int64 = 1
    last_allotment::Int64 = 1
    observations_total::Int64 = 0
    adaptation_updates_total::Int64 = 0
    last_elapsed_ns::Int64 = 0
    elapsed_ns_total::Int64 = 0
    threaded_elapsed_ns_total::Int64 = 0
    serial_elapsed_ns_total::Int64 = 0
    last_classification::Symbol = :none
    last_utilization::Float64 = 1.0
    quantum_length::Int64 = 0
    trim_quanta_budget::Int64 = 0
    quantums_total::Int64 = 0
    quantums_inefficient::Int64 = 0
    quantums_efficient_satisfied::Int64 = 0
    quantums_efficient_deprived::Int64 = 0
    quantums_accounted_proxy::Int64 = 0
    quantums_deductible_proxy::Int64 = 0
end

Base.@kwdef mutable struct PolicyContext
    telemetry::PolicyTelemetry = PolicyTelemetry()
    adaptive_state::Dict{Symbol, AdaptiveControllerState} = Dict{Symbol, AdaptiveControllerState}()
end

const _policy_telemetry_lock = ReentrantLock()
const _policy_context_tls_key = :spaceagora_parallel_policy_context
const _global_policy_context = Ref{PolicyContext}(PolicyContext())

@inline function _active_policy_context()::PolicyContext
    ctx = try
        Base.task_local_storage(_policy_context_tls_key)
    catch
        nothing
    end
    if ctx isa PolicyContext
        return ctx
    end
    return _global_policy_context[]
end

function with_policy_context(f::Function)
    return Base.task_local_storage(_policy_context_tls_key, PolicyContext()) do
        f()
    end
end

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

@inline function _telemetry_bucket(source::Symbol)::Symbol
    if source == :density_callback
        return :density
    elseif source == :control_callback
        return :control
    elseif source == :multibody
        return :multibody
    end
    return :other
end

@inline adaptive_policy_enabled()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_ADAPTIVE", false)
@inline adaptive_window_size()::Int = parse_thread_threshold_env("SPACEAGORA_PARALLEL_POLICY_WINDOW", 8)
@inline adaptive_trim_quanta_budget()::Int = parse_nonnegative_int_env("SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA", 0)
@inline adaptive_bootstrap_threads()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS", true)
@inline adaptive_control_tail_guard()::Bool = parse_bool_env("SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD", false)

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

@inline function _adaptive_state_for(source::Symbol)::AdaptiveControllerState
    ctx = _active_policy_context()
    return get!(ctx.adaptive_state, source) do
        AdaptiveControllerState()
    end
end

function _record_policy_decision!(
    source::Symbol,
    mode::Symbol,
    threshold::Int,
    num_items::Int,
    budget::Int,
    adaptive_enabled::Bool,
    desire::Int,
    allotment::Int,
    outer_active::Bool,
    allow_with_outer::Bool,
    heavy_only::Bool,
    heavy_work::Bool,
    use_threads::Bool
)
    lock(_policy_telemetry_lock) do
        t = _active_policy_context().telemetry
        t.decisions_total += 1
        t.threads_enabled_total += use_threads ? 1 : 0
        t.adaptive_decisions_total += adaptive_enabled ? 1 : 0

        bucket = _telemetry_bucket(source)
        if bucket == :density
            t.density_decisions += 1
            t.density_threads_enabled += use_threads ? 1 : 0
        elseif bucket == :control
            t.control_decisions += 1
            t.control_threads_enabled += use_threads ? 1 : 0
        elseif bucket == :multibody
            t.multibody_decisions += 1
            t.multibody_threads_enabled += use_threads ? 1 : 0
        else
            t.other_decisions += 1
            t.other_threads_enabled += use_threads ? 1 : 0
        end

        t.last_source = source
        t.last_mode = mode
        t.last_threshold = max(1, threshold)
        t.last_num_items = max(0, num_items)
        t.last_budget = max(1, budget)
        t.last_adaptive_enabled = adaptive_enabled
        t.last_desire = max(1, desire)
        t.last_allotment = max(1, allotment)
        t.last_outer_active = outer_active
        t.last_allow_with_outer = allow_with_outer
        t.last_heavy_only = heavy_only
        t.last_heavy_work = heavy_work
        t.last_use_threads = use_threads
    end
    return nothing
end

function reset_policy_telemetry!()
    lock(_policy_telemetry_lock) do
        ctx = _active_policy_context()
        ctx.telemetry = PolicyTelemetry()
        empty!(ctx.adaptive_state)
    end
    return nothing
end

function policy_telemetry_snapshot()
    lock(_policy_telemetry_lock) do
        t = _active_policy_context().telemetry
        quantums_effective = max(0, t.quantums_total - min(t.quantums_total, t.trim_quanta_budget))
        trimmed_accounted = min(t.quantums_accounted_proxy, quantums_effective)
        accounted_fraction_proxy = t.quantums_total == 0 ? 0.0 : t.quantums_accounted_proxy / t.quantums_total
        trimmed_accounted_fraction_proxy = quantums_effective == 0 ? 0.0 : trimmed_accounted / quantums_effective
        return (
            decisions_total=t.decisions_total,
            threads_enabled_total=t.threads_enabled_total,
            adaptive_decisions_total=t.adaptive_decisions_total,
            density_decisions=t.density_decisions,
            density_threads_enabled=t.density_threads_enabled,
            control_decisions=t.control_decisions,
            control_threads_enabled=t.control_threads_enabled,
            multibody_decisions=t.multibody_decisions,
            multibody_threads_enabled=t.multibody_threads_enabled,
            other_decisions=t.other_decisions,
            other_threads_enabled=t.other_threads_enabled,
            last_source=String(t.last_source),
            last_mode=String(t.last_mode),
            last_threshold=t.last_threshold,
            last_num_items=t.last_num_items,
            last_budget=t.last_budget,
            last_adaptive_enabled=t.last_adaptive_enabled,
            last_desire=t.last_desire,
            last_allotment=t.last_allotment,
            last_outer_active=t.last_outer_active,
            last_allow_with_outer=t.last_allow_with_outer,
            last_heavy_only=t.last_heavy_only,
            last_heavy_work=t.last_heavy_work,
            last_use_threads=t.last_use_threads,
            observations_total=t.observations_total,
            adaptation_updates_total=t.adaptation_updates_total,
            last_elapsed_ns=t.last_elapsed_ns,
            elapsed_ns_total=t.elapsed_ns_total,
            threaded_elapsed_ns_total=t.threaded_elapsed_ns_total,
            serial_elapsed_ns_total=t.serial_elapsed_ns_total,
            last_classification=String(t.last_classification),
            last_utilization=t.last_utilization,
            quantum_length=t.quantum_length,
            trim_quanta_budget=t.trim_quanta_budget,
            quantums_total=t.quantums_total,
            quantums_inefficient=t.quantums_inefficient,
            quantums_efficient_satisfied=t.quantums_efficient_satisfied,
            quantums_efficient_deprived=t.quantums_efficient_deprived,
            quantums_accounted_proxy=t.quantums_accounted_proxy,
            quantums_deductible_proxy=t.quantums_deductible_proxy,
            accounted_fraction_proxy=accounted_fraction_proxy,
            trimmed_accounted_fraction_proxy=trimmed_accounted_fraction_proxy
        )
    end
end

@inline function thread_policy_decision(
    num_items::Int;
    mode::Symbol,
    threshold::Int,
    heavy_work::Bool=true,
    heavy_only::Bool=false,
    outer_active::Bool=false,
    allow_with_outer::Bool=true,
    source::Symbol=:other
)
    budget = effective_inner_thread_budget()
    adaptive_enabled = (mode == :auto) && adaptive_policy_enabled()
    bootstrap_threads = adaptive_enabled && adaptive_bootstrap_threads()
    control_tail_guard = adaptive_enabled && adaptive_control_tail_guard()
    desire = 1
    allotment = 1
    if adaptive_enabled
        ρ = adaptive_rho()
        desire_cap = _adaptive_desire_cap(budget, ρ)
        lock(_policy_telemetry_lock) do
            st = _adaptive_state_for(source)
            st.desire = min(max(1, st.desire), desire_cap)
            # Tail guard: avoid a cold-start serial window on obviously parallel workloads.
            if bootstrap_threads && st.desire == 1 && budget > 1 && num_items >= max(1, threshold)
                st.desire = min(desire_cap, 2)
            end
            if control_tail_guard && source == :control_callback && budget > 1 && num_items >= max(1, threshold)
                stable_desire = min(desire_cap, min(budget, max(2, num_items)))
                st.desire = max(st.desire, stable_desire)
            end
            desire = st.desire
        end
        allotment = max(1, min(desire, budget))
    else
        if mode == :on
            desire = max(2, budget)
        elseif mode == :auto
            desire = max(1, budget)
        end
        allotment = if mode == :off
            1
        else
            max(1, min(max(1, desire), budget))
        end
    end

    use_threads =
        if budget <= 1 || num_items <= 1
            false
        elseif mode == :off
            false
        elseif mode == :on
            true
        elseif outer_active && !allow_with_outer
            false
        elseif heavy_only && !heavy_work
            false
        elseif adaptive_enabled
            desire > 1 && num_items >= max(1, threshold)
        else
            num_items >= max(1, threshold)
        end

    allotted = use_threads ? allotment : 1
    _record_policy_decision!(
        source,
        mode,
        threshold,
        num_items,
        budget,
        adaptive_enabled,
        desire,
        allotted,
        outer_active,
        allow_with_outer,
        heavy_only,
        heavy_work,
        use_threads
    )
    return (
        use_threads=use_threads,
        allotment=allotted,
        budget=budget,
        mode=mode,
        threshold=max(1, threshold),
        num_items=max(0, num_items),
        adaptive_enabled=adaptive_enabled,
        desire=max(1, desire)
    )
end

@inline function use_threads_policy(
    num_items::Int;
    mode::Symbol,
    threshold::Int,
    heavy_work::Bool=true,
    heavy_only::Bool=false,
    outer_active::Bool=false,
    allow_with_outer::Bool=true,
    source::Symbol=:other
)::Bool
    return thread_policy_decision(
        num_items;
        mode=mode,
        threshold=threshold,
        heavy_work=heavy_work,
        heavy_only=heavy_only,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=source
    ).use_threads
end

function threaded_foreach(num_items::Int, allotment::Int, f::F) where {F <: Function}
    num_items <= 0 && return nothing
    budget = effective_inner_thread_budget()
    workers = min(num_items, max(1, allotment), budget)
    if workers <= 1 || Threads.nthreads() <= 1
        @inbounds for idx in 1:num_items
            f(idx)
        end
        return nothing
    end
    Threads.@sync for worker_id in 1:workers
        Threads.@spawn begin
            @inbounds for idx in worker_id:workers:num_items
                f(idx)
            end
        end
    end
    return nothing
end

function threaded_foreach(f::F, num_items::Int, allotment::Int) where {F <: Function}
    return threaded_foreach(num_items, allotment, f)
end

function record_policy_observation!(
    source::Symbol;
    mode::Symbol,
    num_items::Int,
    use_threads::Bool,
    elapsed_ns::Integer
)
    budget = effective_inner_thread_budget()
    elapsed_ns_i64 = try
        Int64(elapsed_ns)
    catch
        typemax(Int64)
    end
    elapsed_ns_clamped = max(Int64(0), elapsed_ns_i64)
    adaptive_enabled = (mode == :auto) && adaptive_policy_enabled()

    lock(_policy_telemetry_lock) do
        t = _active_policy_context().telemetry
        t.observations_total += 1
        t.last_elapsed_ns = elapsed_ns_clamped
        t.elapsed_ns_total += elapsed_ns_clamped
        if use_threads
            t.threaded_elapsed_ns_total += elapsed_ns_clamped
        else
            t.serial_elapsed_ns_total += elapsed_ns_clamped
        end

        if !adaptive_enabled
            return nothing
        end

        st = _adaptive_state_for(source)
        ρ = adaptive_rho()
        δ = adaptive_delta()
        L = adaptive_window_size()
        trim_quanta = adaptive_trim_quanta_budget()

        st.desire = min(max(1, st.desire), _adaptive_desire_cap(budget, ρ))
        allotment = use_threads ? max(1, min(st.desire, budget)) : 1
        useful = min(max(0, num_items), allotment)

        st.window_calls += 1
        st.window_allotment_sum += allotment
        st.window_useful_sum += useful
        st.window_deprived_calls += (allotment < st.desire) ? 1 : 0

        if st.window_calls < L
            return nothing
        end

        utilization = st.window_useful_sum / max(1, st.window_allotment_sum)
        efficient = utilization >= δ
        if !efficient
            st.last_classification = :inefficient
            st.desire = max(1, floor(Int, st.desire / ρ))
        elseif st.window_deprived_calls == 0
            st.last_classification = :efficient_satisfied
            st.desire = min(_adaptive_desire_cap(budget, ρ), max(1, ceil(Int, st.desire * ρ)))
        else
            st.last_classification = :efficient_deprived
        end

        st.last_utilization = utilization
        st.window_calls = 0
        st.window_allotment_sum = 0
        st.window_useful_sum = 0.0
        st.window_deprived_calls = 0

        t.adaptation_updates_total += 1
        t.last_classification = st.last_classification
        t.last_utilization = st.last_utilization
        t.last_desire = st.desire
        t.quantum_length = L
        t.trim_quanta_budget = trim_quanta
        t.quantums_total += 1
        if st.last_classification == :inefficient
            t.quantums_inefficient += 1
            t.quantums_deductible_proxy += 1
        elseif st.last_classification == :efficient_satisfied
            t.quantums_efficient_satisfied += 1
            t.quantums_deductible_proxy += 1
        elseif st.last_classification == :efficient_deprived
            t.quantums_efficient_deprived += 1
            t.quantums_accounted_proxy += 1
        else
            t.quantums_deductible_proxy += 1
        end
    end
    return nothing
end

end # module ParallelPolicy
