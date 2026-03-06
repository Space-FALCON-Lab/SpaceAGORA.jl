module ParallelPolicy

using Base.Threads
using TOML

export parse_bool_env, parse_parallel_mode_env, parse_thread_threshold_env
export outer_parallel_active, effective_inner_thread_budget, use_threads_policy
export thread_policy_decision, threaded_foreach, threaded_reduce, threaded_foreach_persistent, with_policy_context
export reset_policy_telemetry!, policy_telemetry_snapshot, record_policy_observation!
export reset_persistent_hint_state!, persistent_hints_state_reset_requested
export hint_layer_stats_snapshot

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
    persistent_hints_enabled::Bool = false
    persistent_hints_loaded::Bool = false
    persistent_hints_updates::Int64 = 0
    persistent_hints_entries::Int64 = 0
    persistent_hints_path::String = ""
    last_signature::String = ""
    last_hint_allotment::Int64 = 1
    last_hint_confidence::Float64 = 0.0
    last_hint_regret_ns::Float64 = 0.0
end

Base.@kwdef mutable struct AdaptiveChoiceStats
    samples::Int64 = 0
    successes::Int64 = 0
    failures::Int64 = 0
    elapsed_sum_ns::Float64 = 0.0
    elapsed_sq_sum_ns::Float64 = 0.0
end

Base.@kwdef mutable struct _HintLayerStatsAccumulator
    signatures::Set{String} = Set{String}()
    choice_count::Int64 = 0
    samples_total::Int64 = 0
    successes_total::Int64 = 0
    failures_total::Int64 = 0
    elapsed_sum_ns::Float64 = 0.0
    elapsed_sq_sum_ns::Float64 = 0.0
    confidence_sum::Float64 = 0.0
    regret_sum_ns::Float64 = 0.0
    signature_metric_count::Int64 = 0
end

Base.@kwdef mutable struct PolicyContext
    telemetry::PolicyTelemetry = PolicyTelemetry()
    adaptive_state::Dict{Symbol, AdaptiveControllerState} = Dict{Symbol, AdaptiveControllerState}()
    decision_signature::Dict{Symbol, String} = Dict{Symbol, String}()
    decision_allotment::Dict{Symbol, Int64} = Dict{Symbol, Int64}()
end

const _policy_telemetry_lock = ReentrantLock()
const _policy_context_tls_key = :spaceagora_parallel_policy_context
const _global_policy_context = Ref{PolicyContext}(PolicyContext())
const _persistent_foreach_lock = ReentrantLock()

Base.@kwdef mutable struct _PersistentHintState
    loaded::Bool = false
    dirty::Bool = false
    path::String = ""
    history::Dict{String, Dict{Int64, AdaptiveChoiceStats}} = Dict{String, Dict{Int64, AdaptiveChoiceStats}}()
end

const _persistent_hint_lock = ReentrantLock()
const _persistent_hint_state = Ref{_PersistentHintState}(_PersistentHintState())
const _persistent_hint_atexit_registered = Ref(false)

Base.@kwdef mutable struct _PersistentForeachPool
    workers::Int
    request_channels::Vector{Channel{Any}}
    done_channel::Channel{Any}
    run_lock::ReentrantLock = ReentrantLock()
end

const _persistent_foreach_pools = Dict{Tuple{UInt, Symbol}, _PersistentForeachPool}()

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

@inline function _policy_scope_id(ctx::PolicyContext)::UInt
    return UInt(objectid(ctx))
end

@inline function _active_policy_scope_id()::UInt
    return _policy_scope_id(_active_policy_context())
end

function _persistent_foreach_worker_loop(
    worker_id::Int,
    request_channel::Channel{Any},
    done_channel::Channel{Any}
)::Nothing
    while true
        request = take!(request_channel)
        if request === :stop
            return nothing
        end

        captured = nothing
        try
            num_items = request.num_items
            active_workers = request.active_workers
            scheduler = request.scheduler
            chunk = request.chunk
            f = request.f
            if scheduler == :dynamic
                next_index = request.next_index
                @inbounds while true
                    start_idx = Threads.atomic_add!(next_index, chunk)
                    start_idx > num_items && break
                    stop_idx = min(num_items, start_idx + chunk - 1)
                    for idx in start_idx:stop_idx
                        f(idx)
                    end
                end
            else
                @inbounds for idx in worker_id:active_workers:num_items
                    f(idx)
                end
            end
        catch err
            captured = Base.CapturedException(err, catch_backtrace())
        end
        put!(done_channel, captured)
    end
end

function _create_persistent_foreach_pool(workers::Int)::_PersistentForeachPool
    workers = max(2, workers)
    request_channels = [Channel{Any}(1) for _ in 1:workers]
    done_channel = Channel{Any}(workers)
    pool = _PersistentForeachPool(
        workers=workers,
        request_channels=request_channels,
        done_channel=done_channel
    )
    @inbounds for worker_id in 1:workers
        Threads.@spawn _persistent_foreach_worker_loop(
            worker_id,
            request_channels[worker_id],
            done_channel
        )
    end
    return pool
end

function _shutdown_persistent_foreach_pool!(pool::_PersistentForeachPool)::Nothing
    lock(pool.run_lock) do
        @inbounds for channel in pool.request_channels
            put!(channel, :stop)
        end
    end
    return nothing
end

function _destroy_persistent_foreach_scope!(scope_id::UInt)::Nothing
    pools = _PersistentForeachPool[]
    lock(_persistent_foreach_lock) do
        stale_keys = Tuple{UInt, Symbol}[]
        for (key, pool) in _persistent_foreach_pools
            if key[1] == scope_id
                push!(stale_keys, key)
                push!(pools, pool)
            end
        end
        @inbounds for key in stale_keys
            delete!(_persistent_foreach_pools, key)
        end
    end
    @inbounds for pool in pools
        _shutdown_persistent_foreach_pool!(pool)
    end
    return nothing
end

function with_policy_context(f::Function)
    ctx = PolicyContext()
    scope_id = _policy_scope_id(ctx)
    return Base.task_local_storage(_policy_context_tls_key, ctx) do
        try
            return f()
        finally
            _destroy_persistent_foreach_scope!(scope_id)
        end
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

@inline function _hint_entry_count(state::_PersistentHintState)::Int
    entries = 0
    for bucket in values(state.history)
        entries += length(bucket)
    end
    return entries
end

@inline function _hint_stats_payload(stats::AdaptiveChoiceStats)::Dict{String, Any}
    return Dict{String, Any}(
        "samples" => max(0, Int(stats.samples)),
        "successes" => max(0, Int(stats.successes)),
        "failures" => max(0, Int(stats.failures)),
        "elapsed_sum_ns" => max(0.0, Float64(stats.elapsed_sum_ns)),
        "elapsed_sq_sum_ns" => max(0.0, Float64(stats.elapsed_sq_sum_ns))
    )
end

@inline function _hint_payload_stats(payload)::Union{Nothing, AdaptiveChoiceStats}
    payload isa AbstractDict || return nothing
    samples = try
        max(0, Int(get(payload, "samples", 0)))
    catch
        0
    end
    successes = try
        max(0, Int(get(payload, "successes", 0)))
    catch
        0
    end
    failures = try
        max(0, Int(get(payload, "failures", 0)))
    catch
        0
    end
    elapsed_sum_ns = try
        max(0.0, Float64(get(payload, "elapsed_sum_ns", 0.0)))
    catch
        0.0
    end
    elapsed_sq_sum_ns = try
        max(0.0, Float64(get(payload, "elapsed_sq_sum_ns", 0.0)))
    catch
        0.0
    end
    samples <= 0 && return nothing
    successes = min(samples, successes)
    failures = min(samples - successes, failures)
    return AdaptiveChoiceStats(
        samples=samples,
        successes=successes,
        failures=failures,
        elapsed_sum_ns=elapsed_sum_ns,
        elapsed_sq_sum_ns=elapsed_sq_sum_ns
    )
end

function _load_persistent_hint_state_locked!()::Nothing
    state = _persistent_hint_state[]
    if state.loaded
        return nothing
    end
    state.loaded = true
    state.path = _persistent_hint_path()
    if persistent_hints_state_reset_requested()
        # Cold-start mode: intentionally skip loading prior persisted hints.
        state.dirty = false
        empty!(state.history)
        return nothing
    end
    if !persistent_hints_enabled()
        return nothing
    end
    if !isfile(state.path)
        return nothing
    end
    parsed = try
        TOML.parsefile(state.path)
    catch
        return nothing
    end
    rows = get(parsed, "history", Any[])
    rows isa AbstractVector || return nothing
    for row in rows
        row isa AbstractDict || continue
        signature = strip(String(get(row, "signature", "")))
        isempty(signature) && continue
        allotment = try
            Int64(get(row, "allotment", 0))
        catch
            0
        end
        allotment > 0 || continue
        stats = _hint_payload_stats(get(row, "stats", nothing))
        stats === nothing && continue
        bucket = get!(state.history, signature) do
            Dict{Int64, AdaptiveChoiceStats}()
        end
        existing = get!(bucket, allotment) do
            AdaptiveChoiceStats()
        end
        existing.samples += stats.samples
        existing.successes += stats.successes
        existing.failures += stats.failures
        existing.elapsed_sum_ns += stats.elapsed_sum_ns
        existing.elapsed_sq_sum_ns += stats.elapsed_sq_sum_ns
    end
    return nothing
end

function _save_persistent_hint_state_locked!()::Nothing
    state = _persistent_hint_state[]
    if !state.loaded || !state.dirty || !persistent_hints_persist_enabled()
        return nothing
    end
    rows = Dict{String, Any}[]
    for signature in sort!(collect(keys(state.history)))
        bucket = state.history[signature]
        isempty(bucket) && continue
        for allotment in sort!(collect(keys(bucket)))
            stats = bucket[allotment]
            stats.samples > 0 || continue
            push!(rows, Dict{String, Any}(
                "signature" => signature,
                "allotment" => Int(allotment),
                "stats" => _hint_stats_payload(stats)
            ))
        end
    end
    payload = Dict{String, Any}(
        "schema_version" => 1,
        "history" => rows
    )
    mkpath(dirname(state.path))
    tmp_path = state.path * ".tmp"
    open(tmp_path, "w") do io
        TOML.print(io, payload)
    end
    mv(tmp_path, state.path; force=true)
    state.dirty = false
    return nothing
end

function _ensure_persistent_hint_state_loaded!()::Nothing
    lock(_persistent_hint_lock) do
        _load_persistent_hint_state_locked!()
        if persistent_hints_enabled() && !_persistent_hint_atexit_registered[]
            _persistent_hint_atexit_registered[] = true
            atexit(() -> begin
                lock(_persistent_hint_lock) do
                    _save_persistent_hint_state_locked!()
                end
                nothing
            end)
        end
    end
    return nothing
end

function reset_persistent_hint_state!()::Nothing
    lock(_persistent_hint_lock) do
        state = _persistent_hint_state[]
        state.loaded = false
        state.dirty = false
        state.path = ""
        empty!(state.history)
    end
    return nothing
end

@inline function _hint_bucket(v::Int)::String
    if v <= 1
        return "1"
    elseif v <= 2
        return "2"
    elseif v <= 4
        return "3_4"
    elseif v <= 8
        return "5_8"
    elseif v <= 16
        return "9_16"
    end
    return "17p"
end

@inline function _hint_workload_signature(
    source::Symbol,
    num_items::Int,
    threshold::Int,
    budget::Int,
    outer_active::Bool,
    heavy_only::Bool,
    heavy_work::Bool
)::String
    profile = _safe_token(get(ENV, "SPACEAGORA_PARALLEL_PROFILE", "default"))
    machine = _safe_token(get(ENV, "SPACEAGORA_PERF_MACHINE_LABEL", "default"))
    return join((
        "profile=$profile",
        "machine=$machine",
        "src=$(String(source))",
        "items=$(_hint_bucket(max(0, num_items)))",
        "thr=$(_hint_bucket(max(1, threshold)))",
        "budget=$(_hint_bucket(max(1, budget)))",
        "outer=$(outer_active ? "1" : "0")",
        "heavy_only=$(heavy_only ? "1" : "0")",
        "heavy=$(heavy_work ? "1" : "0")"
    ), "|")
end

@inline function _hint_candidate_allotments(num_items::Int, budget::Int)::Vector{Int64}
    cap = max(1, min(max(1, num_items), max(1, budget)))
    candidates = Int64[1]
    if cap >= 2
        push!(candidates, 2)
    end
    if cap >= 4
        push!(candidates, 4)
    end
    if cap >= 8
        push!(candidates, 8)
    end
    push!(candidates, Int64(cld(cap, 2)))
    push!(candidates, Int64(cap))
    sort!(unique!(candidates))
    return candidates
end

@inline function _hint_mean_and_width(
    stats::AdaptiveChoiceStats,
    total_samples::Int64,
    explore_c::Float64
)::Tuple{Float64, Float64}
    n = max(Int64(1), stats.samples)
    mean_ns = stats.elapsed_sum_ns / n
    width = explore_c * sqrt(log(max(2.0, Float64(total_samples))) / n)
    return mean_ns, width
end

@inline function _hint_samples(bucket, candidate::Int64)::Int64
    stats = get(bucket, candidate, nothing)
    if stats isa AdaptiveChoiceStats
        return max(Int64(0), stats.samples)
    end
    return Int64(0)
end

function _hint_choose_allotment(
    signature::String,
    candidates::Vector{Int64}
)::NamedTuple{(:allotment, :confidence, :regret_ns, :samples, :exploring), Tuple{Int64, Float64, Float64, Int64, Bool}}
    _ensure_persistent_hint_state_loaded!()
    candidate_pool = isempty(candidates) ? Int64[1] : sort!(unique!(Int64[max(Int64(1), c) for c in candidates]))
    if !persistent_hints_enabled()
        return (allotment=Int64(1), confidence=0.0, regret_ns=0.0, samples=Int64(0), exploring=false)
    end
    lock(_persistent_hint_lock) do
        state = _persistent_hint_state[]
        bucket = get(state.history, signature, nothing)
        if bucket === nothing || isempty(bucket)
            return (
                allotment=first(candidate_pool),
                confidence=Inf,
                regret_ns=0.0,
                samples=Int64(0),
                exploring=true
            )
        end
        total_samples = Int64(0)
        for c in candidate_pool
            stats = get(bucket, c, nothing)
            stats isa AdaptiveChoiceStats || continue
            total_samples += stats.samples
        end
        total_samples = max(Int64(0), total_samples)

        explore_c = persistent_hints_exploration()
        min_samples = Int64(persistent_hints_min_samples())
        explore_candidate = Int64(0)
        explore_samples = typemax(Int64)
        for c in candidate_pool
            samples = _hint_samples(bucket, c)
            if samples < min_samples && (samples < explore_samples || (samples == explore_samples && (explore_candidate == 0 || c < explore_candidate)))
                explore_candidate = c
                explore_samples = samples
            end
        end
        if explore_candidate > 0
            stats = get(bucket, explore_candidate, nothing)
            confidence = if stats isa AdaptiveChoiceStats && stats.samples > 0 && total_samples > 0
                _, width = _hint_mean_and_width(stats, total_samples, explore_c)
                width
            else
                Inf
            end
            return (
                allotment=explore_candidate,
                confidence=confidence,
                regret_ns=0.0,
                samples=max(Int64(0), explore_samples),
                exploring=true
            )
        end

        best_allotment = Int64(1)
        best_score = Inf
        best_width = 0.0
        best_mean = Inf
        best_samples = Int64(0)
        known_means = Float64[]
        for c in candidate_pool
            stats = get(bucket, c, nothing)
            stats isa AdaptiveChoiceStats || continue
            stats.samples > 0 || continue
            mean_ns, width = _hint_mean_and_width(stats, total_samples, explore_c)
            score = mean_ns - width
            push!(known_means, mean_ns)
            if score < best_score
                best_score = score
                best_allotment = c
                best_width = width
                best_mean = mean_ns
                best_samples = max(Int64(0), stats.samples)
            end
        end
        if isempty(known_means)
            fallback = first(candidate_pool)
            return (allotment=fallback, confidence=Inf, regret_ns=0.0, samples=Int64(0), exploring=true)
        end
        best_observed = minimum(known_means)
        regret = max(0.0, best_mean - best_observed)
        return (
            allotment=best_allotment,
            confidence=best_width,
            regret_ns=regret,
            samples=best_samples,
            exploring=false
        )
    end
end

function _hint_record_observation!(
    signature::String,
    allotment::Int64,
    elapsed_ns::Int64;
    success::Bool
)::Nothing
    _ensure_persistent_hint_state_loaded!()
    if !persistent_hints_enabled()
        return nothing
    end
    allotment > 0 || return nothing
    elapsed = max(0.0, Float64(elapsed_ns))
    lock(_persistent_hint_lock) do
        state = _persistent_hint_state[]
        bucket = get!(state.history, signature) do
            Dict{Int64, AdaptiveChoiceStats}()
        end
        stats = get!(bucket, allotment) do
            AdaptiveChoiceStats()
        end
        stats.samples += 1
        stats.successes += success ? 1 : 0
        stats.failures += success ? 0 : 1
        stats.elapsed_sum_ns += elapsed
        stats.elapsed_sq_sum_ns += elapsed^2
        state.dirty = true
    end
    return nothing
end

@inline function _hint_signature_value(signature::String, key::String)::Union{Nothing, String}
    prefix = string(key, "=")
    n = lastindex(prefix)
    for token in split(signature, '|')
        startswith(token, prefix) || continue
        if n >= lastindex(token)
            return ""
        end
        return String(token[(n + 1):end])
    end
    return nothing
end

function hint_layer_stats_snapshot(;
    profile::Union{Nothing, AbstractString}=nothing,
    machine::Union{Nothing, AbstractString}=nothing
)
    _ensure_persistent_hint_state_loaded!()
    profile_filter = profile === nothing ? nothing : _safe_token(String(profile))
    machine_filter = machine === nothing ? nothing : _safe_token(String(machine))
    explore_c = persistent_hints_exploration()
    min_samples = Int64(persistent_hints_min_samples())
    rows = NamedTuple[]

    lock(_persistent_hint_lock) do
        state = _persistent_hint_state[]
        if !state.loaded || isempty(state.history)
            return nothing
        end

        acc = Dict{Tuple{String, String, Symbol}, _HintLayerStatsAccumulator}()
        for (signature, bucket) in state.history
            bucket isa Dict{Int64, AdaptiveChoiceStats} || continue
            profile_token = something(_hint_signature_value(signature, "profile"), "unknown")
            machine_token = something(_hint_signature_value(signature, "machine"), "unknown")
            source_token = something(_hint_signature_value(signature, "src"), "other")
            source_sym = Symbol(source_token)
            if !(profile_filter === nothing || profile_token == profile_filter)
                continue
            end
            if !(machine_filter === nothing || machine_token == machine_filter)
                continue
            end

            entries = Tuple{Int64, AdaptiveChoiceStats}[]
            total_signature_samples = Int64(0)
            for (allotment, stats) in bucket
                stats.samples > 0 || continue
                push!(entries, (max(Int64(1), allotment), stats))
                total_signature_samples += stats.samples
            end
            isempty(entries) && continue

            key = (profile_token, machine_token, source_sym)
            layer = get!(acc, key) do
                _HintLayerStatsAccumulator()
            end
            push!(layer.signatures, signature)

            best_score = Inf
            best_mean = Inf
            best_width = 0.0
            observed_best_mean = Inf
            for (_, stats) in entries
                layer.choice_count += 1
                layer.samples_total += stats.samples
                layer.successes_total += stats.successes
                layer.failures_total += stats.failures
                layer.elapsed_sum_ns += stats.elapsed_sum_ns
                layer.elapsed_sq_sum_ns += stats.elapsed_sq_sum_ns
                mean_ns, width = _hint_mean_and_width(
                    stats,
                    max(Int64(1), total_signature_samples),
                    explore_c
                )
                observed_best_mean = min(observed_best_mean, mean_ns)
                score = mean_ns - width
                if score < best_score
                    best_score = score
                    best_mean = mean_ns
                    best_width = width
                end
            end
            if isfinite(best_mean) && isfinite(observed_best_mean)
                layer.regret_sum_ns += max(0.0, best_mean - observed_best_mean)
                layer.confidence_sum += max(0.0, best_width)
                layer.signature_metric_count += 1
            end
        end

        ordered_keys = sort(collect(keys(acc)); by=k -> (k[1], k[2], String(k[3])))
        for key in ordered_keys
            layer = acc[key]
            n = max(Int64(1), layer.samples_total)
            elapsed_mean_ns = layer.elapsed_sum_ns / n
            elapsed_var_ns = max(0.0, layer.elapsed_sq_sum_ns / n - elapsed_mean_ns^2)
            metric_n = max(Int64(1), layer.signature_metric_count)
            push!(rows, (
                profile=key[1],
                machine=key[2],
                layer=String(key[3]),
                signature_count=Int64(length(layer.signatures)),
                choice_count=layer.choice_count,
                samples_total=layer.samples_total,
                successes_total=layer.successes_total,
                failures_total=layer.failures_total,
                elapsed_mean_ns=elapsed_mean_ns,
                elapsed_std_ns=sqrt(elapsed_var_ns),
                confidence_mean=layer.confidence_sum / metric_n,
                regret_mean_ns=layer.regret_sum_ns / metric_n,
                exploration_c=explore_c,
                min_samples=min_samples,
                state_path=state.path
            ))
        end
        return nothing
    end

    return rows
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
    use_threads::Bool,
    signature::String,
    hint_allotment::Int64,
    hint_confidence::Float64,
    hint_regret_ns::Float64,
    hints_loaded::Bool,
    hints_entries::Int64
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
        t.persistent_hints_enabled = persistent_hints_enabled()
        t.persistent_hints_loaded = hints_loaded
        t.persistent_hints_entries = hints_entries
        t.persistent_hints_path = _persistent_hint_state[].path
        t.last_signature = signature
        t.last_hint_allotment = max(1, hint_allotment)
        t.last_hint_confidence = max(0.0, hint_confidence)
        t.last_hint_regret_ns = max(0.0, hint_regret_ns)
    end
    return nothing
end

function reset_policy_telemetry!()
    lock(_policy_telemetry_lock) do
        ctx = _active_policy_context()
        ctx.telemetry = PolicyTelemetry()
        empty!(ctx.adaptive_state)
        empty!(ctx.decision_signature)
        empty!(ctx.decision_allotment)
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
            persistent_hints_enabled=t.persistent_hints_enabled,
            persistent_hints_loaded=t.persistent_hints_loaded,
            persistent_hints_updates=t.persistent_hints_updates,
            persistent_hints_entries=t.persistent_hints_entries,
            persistent_hints_path=t.persistent_hints_path,
            last_signature=t.last_signature,
            last_hint_allotment=t.last_hint_allotment,
            last_hint_confidence=t.last_hint_confidence,
            last_hint_regret_ns=t.last_hint_regret_ns,
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
    measured_reward = adaptive_enabled && persistent_hints_enabled() && adaptive_measured_reward_enabled()
    bootstrap_threads = adaptive_enabled && adaptive_bootstrap_threads()
    control_tail_guard = adaptive_enabled && adaptive_control_tail_guard()
    signature = ""
    hint_allotment = Int64(1)
    hint_confidence = 0.0
    hint_regret_ns = 0.0
    hints_loaded = false
    hints_entries = 0
    desire = 1
    allotment = 1
    if adaptive_enabled
        signature = _hint_workload_signature(
            source,
            num_items,
            threshold,
            budget,
            outer_active,
            heavy_only,
            heavy_work
        )
        hint = _hint_choose_allotment(signature, _hint_candidate_allotments(num_items, budget))
        hint_allotment = hint.allotment
        hint_confidence = hint.confidence
        hint_regret_ns = hint.regret_ns
        lock(_persistent_hint_lock) do
            hints_loaded = _persistent_hint_state[].loaded
            hints_entries = _hint_entry_count(_persistent_hint_state[])
        end
        ρ = adaptive_rho()
        desire_cap = _adaptive_desire_cap(budget, ρ)
        lock(_policy_telemetry_lock) do
            st = _adaptive_state_for(source)
            st.desire = min(max(1, st.desire), desire_cap)
            if measured_reward
                # Measured-reward mode: drive desire directly from elapsed-time hint choice.
                st.desire = min(desire_cap, max(1, Int(hint_allotment)))
            else
                # Tail guard: avoid a cold-start serial window on obviously parallel workloads.
                if bootstrap_threads && st.desire == 1 && budget > 1 && num_items >= max(1, threshold)
                    st.desire = min(desire_cap, 2)
                end
                if control_tail_guard && source == :control_callback && budget > 1 && num_items >= max(1, threshold)
                    stable_desire = min(desire_cap, min(budget, max(2, num_items)))
                    st.desire = max(st.desire, stable_desire)
                end
                if hint_allotment > 1
                    # Blend persisted hint with live AIMD state; this reuses past wins without hard pinning.
                    blended = max(st.desire, Int(hint_allotment))
                    st.desire = min(desire_cap, max(1, blended))
                end
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
        use_threads,
        signature,
        hint_allotment,
        hint_confidence,
        hint_regret_ns,
        hints_loaded,
        hints_entries
    )
    lock(_policy_telemetry_lock) do
        ctx = _active_policy_context()
        ctx.decision_signature[source] = signature
        ctx.decision_allotment[source] = Int64(allotted)
    end
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
    workers = _thread_worker_count(num_items, allotment)
    if workers <= 1 || Threads.nthreads() <= 1
        @inbounds for idx in 1:num_items
            f(idx)
        end
        return nothing
    end
    scheduler = inner_scheduler_mode()
    if scheduler == :dynamic
        chunk = inner_dynamic_chunk_size()
        next_index = Threads.Atomic{Int}(1)
        Threads.@sync for _ in 1:workers
            Threads.@spawn begin
                @inbounds while true
                    start_idx = Threads.atomic_add!(next_index, chunk)
                    start_idx > num_items && break
                    stop_idx = min(num_items, start_idx + chunk - 1)
                    for idx in start_idx:stop_idx
                        f(idx)
                    end
                end
            end
        end
    else
        Threads.@sync for worker_id in 1:workers
            Threads.@spawn begin
                @inbounds for idx in worker_id:workers:num_items
                    f(idx)
                end
            end
        end
    end
    return nothing
end

function threaded_foreach(f::F, num_items::Int, allotment::Int) where {F <: Function}
    return threaded_foreach(num_items, allotment, f)
end

@inline callback_persistent_workers_enabled()::Bool = parse_bool_env("SPACEAGORA_CALLBACK_PERSISTENT_WORKERS", true)

@inline function _persistent_pool_key(source::Symbol)::Tuple{UInt, Symbol}
    return (_active_policy_scope_id(), source)
end

function _persistent_pool_for(source::Symbol)::_PersistentForeachPool
    key = _persistent_pool_key(source)
    lock(_persistent_foreach_lock) do
        return get!(_persistent_foreach_pools, key) do
            _create_persistent_foreach_pool(_default_thread_pool_size())
        end
    end
end

function _threaded_foreach_persistent!(
    pool::_PersistentForeachPool,
    num_items::Int,
    workers::Int,
    f::F
) where {F <: Function}
    scheduler = inner_scheduler_mode()
    chunk = inner_dynamic_chunk_size()
    next_index = Threads.Atomic{Int}(1)
    lock(pool.run_lock) do
        @inbounds for worker_id in 1:workers
            put!(
                pool.request_channels[worker_id],
                (
                    num_items=num_items,
                    active_workers=workers,
                    scheduler=scheduler,
                    chunk=chunk,
                    next_index=next_index,
                    f=f
                )
            )
        end
        first_error = nothing
        @inbounds for _ in 1:workers
            captured = take!(pool.done_channel)
            if !(captured === nothing) && first_error === nothing
                first_error = captured
            end
        end
        if !(first_error === nothing)
            throw(first_error)
        end
    end
    return nothing
end

function threaded_foreach_persistent(
    source::Symbol,
    num_items::Int,
    allotment::Int,
    f::F
) where {F <: Function}
    num_items <= 0 && return nothing
    workers = _thread_worker_count(num_items, allotment)
    if workers <= 1 || !callback_persistent_workers_enabled()
        return threaded_foreach(num_items, allotment, f)
    end
    pool = _persistent_pool_for(source)
    return _threaded_foreach_persistent!(pool, num_items, workers, f)
end

function threaded_foreach_persistent(
    f::F,
    source::Symbol,
    num_items::Int,
    allotment::Int
) where {F <: Function}
    return threaded_foreach_persistent(source, num_items, allotment, f)
end

@inline function _thread_worker_count(num_items::Int, allotment::Int)::Int
    num_items <= 0 && return 1
    budget = effective_inner_thread_budget()
    workers = min(num_items, max(1, allotment), budget)
    if workers <= 1 || Threads.nthreads() <= 1
        return 1
    end
    return workers
end

@inline function thread_worker_count(num_items::Int, allotment::Int)::Int
    return _thread_worker_count(num_items, allotment)
end

function threaded_foreach_worker(num_items::Int, allotment::Int, f::F) where {F <: Function}
    num_items <= 0 && return nothing
    workers = _thread_worker_count(num_items, allotment)
    if workers <= 1
        @inbounds for idx in 1:num_items
            f(1, idx)
        end
        return nothing
    end
    scheduler = inner_scheduler_mode()
    if scheduler == :dynamic
        chunk = inner_dynamic_chunk_size()
        next_index = Threads.Atomic{Int}(1)
        Threads.@sync for worker_id in 1:workers
            Threads.@spawn begin
                @inbounds while true
                    start_idx = Threads.atomic_add!(next_index, chunk)
                    start_idx > num_items && break
                    stop_idx = min(num_items, start_idx + chunk - 1)
                    for idx in start_idx:stop_idx
                        f(worker_id, idx)
                    end
                end
            end
        end
    else
        Threads.@sync for worker_id in 1:workers
            Threads.@spawn begin
                @inbounds for idx in worker_id:workers:num_items
                    f(worker_id, idx)
                end
            end
        end
    end
    return nothing
end

function threaded_foreach_worker(f::F, num_items::Int, allotment::Int) where {F <: Function}
    return threaded_foreach_worker(num_items, allotment, f)
end

function threaded_reduce(
    num_items::Int,
    allotment::Int,
    init::I,
    body!::B,
    combine!::C
) where {I <: Function, B <: Function, C <: Function}
    workers = _thread_worker_count(num_items, allotment)
    acc0 = init()
    if num_items <= 0
        return acc0
    end
    if workers <= 1
        @inbounds for idx in 1:num_items
            body!(acc0, idx)
        end
        return acc0
    end

    partials = Vector{typeof(acc0)}(undef, workers)
    partials[1] = acc0
    scheduler = inner_scheduler_mode()
    if scheduler == :dynamic
        chunk = inner_dynamic_chunk_size()
        next_index = Threads.Atomic{Int}(1)
        Threads.@sync for worker_id in 1:workers
            Threads.@spawn begin
                local_acc = worker_id == 1 ? partials[1] : init()
                @inbounds while true
                    start_idx = Threads.atomic_add!(next_index, chunk)
                    start_idx > num_items && break
                    stop_idx = min(num_items, start_idx + chunk - 1)
                    for idx in start_idx:stop_idx
                        body!(local_acc, idx)
                    end
                end
                partials[worker_id] = local_acc
            end
        end
    else
        Threads.@sync for worker_id in 1:workers
            Threads.@spawn begin
                local_acc = worker_id == 1 ? partials[1] : init()
                @inbounds for idx in worker_id:workers:num_items
                    body!(local_acc, idx)
                end
                partials[worker_id] = local_acc
            end
        end
    end

    result = partials[1]
    @inbounds for worker_id in 2:workers
        combine!(result, partials[worker_id])
    end
    return result
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
    measured_reward = adaptive_enabled && persistent_hints_enabled() && adaptive_measured_reward_enabled()
    hint_signature = ""
    hint_allotment = Int64(1)

    lock(_policy_telemetry_lock) do
        ctx = _active_policy_context()
        t = ctx.telemetry
        t.observations_total += 1
        t.last_elapsed_ns = elapsed_ns_clamped
        t.elapsed_ns_total += elapsed_ns_clamped
        if use_threads
            t.threaded_elapsed_ns_total += elapsed_ns_clamped
        else
            t.serial_elapsed_ns_total += elapsed_ns_clamped
        end

        hint_signature = get(ctx.decision_signature, source, "")
        hint_allotment = get(ctx.decision_allotment, source, use_threads ? Int64(max(1, min(budget, num_items))) : Int64(1))
        t.last_signature = hint_signature

        if !adaptive_enabled
            return nothing
        end

        st = _adaptive_state_for(source)
        if measured_reward
            st.desire = max(Int64(1), hint_allotment)
            st.last_classification = :measured_reward
            st.last_utilization = use_threads ? 1.0 : 0.0
            st.window_calls = 0
            st.window_allotment_sum = 0
            st.window_useful_sum = 0.0
            st.window_deprived_calls = 0

            t.adaptation_updates_total += 1
            t.last_classification = st.last_classification
            t.last_utilization = st.last_utilization
            t.last_desire = st.desire
            t.quantum_length = 1
            t.trim_quanta_budget = 0
            t.quantums_total += 1
            t.quantums_accounted_proxy += 1
            return nothing
        end

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

    if adaptive_enabled && persistent_hints_enabled() && !isempty(hint_signature)
        _hint_record_observation!(
            hint_signature,
            max(Int64(1), hint_allotment),
            elapsed_ns_clamped;
            success=(elapsed_ns_clamped > 0)
        )
        hints_loaded = false
        hints_entries = 0
        hints_path = ""
        lock(_persistent_hint_lock) do
            hints_loaded = _persistent_hint_state[].loaded
            hints_entries = _hint_entry_count(_persistent_hint_state[])
            hints_path = _persistent_hint_state[].path
        end
        lock(_policy_telemetry_lock) do
            t = _active_policy_context().telemetry
            t.persistent_hints_updates += 1
            t.persistent_hints_loaded = hints_loaded
            t.persistent_hints_entries = hints_entries
            t.persistent_hints_path = hints_path
        end
    end
    return nothing
end

end # module ParallelPolicy
