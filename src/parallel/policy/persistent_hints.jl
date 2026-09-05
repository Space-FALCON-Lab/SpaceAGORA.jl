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
    explore_c::Float64;
    scaled::Bool=false
)::Tuple{Float64, Float64}
    n = max(Int64(1), stats.samples)
    mean_ns = stats.elapsed_sum_ns / n
    width = explore_c * sqrt(log(max(2.0, Float64(total_samples))) / n)
    if scaled
        # Scale by the arm's own spread, as `_candidate_confidence_width` in the
        # outer-route selector already does. Unscaled, the width is
        # `c * sqrt(log N / n)` in RAW NANOSECONDS -- one or two ns -- against
        # means measured in microseconds to milliseconds, so the bound equalled
        # the mean and the chooser was a plain argmin over two-sample means.
        # `confidence` reported through telemetry was likewise meaningless.
        # Behind SPACEAGORA_PARALLEL_POLICY_V2 so the shipped chooser stays
        # reproducible for the paired comparison.
        mean_sq_ns = stats.elapsed_sq_sum_ns / n
        std_ns = sqrt(max(0.0, mean_sq_ns - mean_ns^2))
        width *= std_ns
    end
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
    candidates::Vector{Int64};
    # Scale the confidence width by each arm's measured spread; see
    # _hint_mean_and_width. Off is the shipped chooser.
    scaled_width::Bool=false
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
            # WIDEST, not narrowest, on a cold signature.
            #
            # `candidate_pool` is sorted ascending, so `first` handed back the
            # narrowest allotment -- typically 1 -- and the exploration loop
            # below then walks the ladder upward. An adaptive profile therefore
            # began every unseen workload effectively serial on the region it
            # was deciding about, and paid the whole ramp again each time
            # re-exploration reset it.
            #
            # The non-adaptive path answers `max(1, budget)` immediately, and on
            # the workload where this was found that answer is right by a factor
            # of 2.8: measured on atmo256_gram_surrogate_10min, threading the
            # density callback takes the RHS from 1338 to 474 us, and the
            # adaptive profile captured only about half of it (954 us) purely
            # because it started narrow and cycled.
            #
            # Starting at the widest candidate makes an unseen workload no worse
            # than the static profiles at the first call, and leaves the layer
            # free to narrow on measured evidence. The risk is inverted rather
            # than removed -- a workload whose optimum is narrow now pays until
            # it learns -- but that cost is bounded by one window, where the ramp
            # cost was paid on every re-exploration.
            return (
                allotment=last(candidate_pool),
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
                _, width = _hint_mean_and_width(stats, total_samples, explore_c; scaled=scaled_width)
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
            mean_ns, width = _hint_mean_and_width(stats, total_samples, explore_c; scaled=scaled_width)
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

# What one hint consultation costs ON THIS MACHINE, measured rather than assumed.
#
# The hint layer is worth consulting when the work a decision guards is large
# compared with the cost of consulting it. That comparison has two sides and
# both are machine-dependent: a lookup costs 183 ns on this repo's 12-core
# reference box and will cost something else elsewhere, and the guarded work
# varies by four orders of magnitude across the workloads in the benchmark
# catalog. Hard-coding either side ports badly, which is why this is probed at
# first use instead.
#
# Measured once per process against a synthetic signature that hits the same
# code path a real consultation does -- candidate enumeration, the store lookup
# and the lock -- and then cached. The probe itself costs roughly a quarter of a
# millisecond, once, against a solve measured in seconds.
const _HINT_OVERHEAD_NS = Ref{Float64}(-1.0)

function hint_overhead_ns()::Float64
    _HINT_OVERHEAD_NS[] >= 0.0 && return _HINT_OVERHEAD_NS[]
    probe_sig = _hint_workload_signature(:_probe, 8, 1, 8, false, false, true)
    candidates = _hint_candidate_allotments(8, 8)
    best = Inf
    try
        for _ in 1:200                      # warm the path
            _hint_choose_allotment(probe_sig, candidates)
        end
        for _ in 1:5                        # min of five blocks of 200
            t0 = time_ns()
            for _ in 1:200
                _hint_choose_allotment(probe_sig, candidates)
            end
            sample = (time_ns() - t0) / 200
            sample < best && (best = sample)
        end
    catch
        # A probe failure must not disable routing. Fall back to a value that
        # keeps the layer enabled for anything but the very cheapest regions.
        best = 200.0
    end
    _HINT_OVERHEAD_NS[] = isfinite(best) ? max(1.0, best) : 200.0
    return _HINT_OVERHEAD_NS[]
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
