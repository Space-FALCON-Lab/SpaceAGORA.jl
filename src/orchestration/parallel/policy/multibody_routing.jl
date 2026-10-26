@inline _parse_bool_env(name::String, default::Bool)::Bool = ParallelPolicy.parse_bool_env(name, default)

# SPACEAGORA_MULTIBODY_PARALLEL, read once per solve. _multibody_thread_decision
# runs per satellite per RHS call from two effector sites, and the ENV read
# (getenv + a fresh String, then strip/lowercase) was the one cost left in it
# after the forced-false short-circuits: one per satellite per call on every
# route. The engine refreshes the cache at each solve start
# (refresh_multibody_parallel_mode!), which is where a harness `withenv` can
# have changed the value; between solves it is a Ref read.
const _MULTIBODY_MODE_CACHE = Ref{Union{Nothing, Symbol}}(nothing)

"""
    refresh_multibody_parallel_mode!() -> Symbol

Re-read `SPACEAGORA_MULTIBODY_PARALLEL` into the per-solve cache the
per-satellite thread decision consults. Called by the engine at solve start.
"""
function refresh_multibody_parallel_mode!()::Symbol
    mode = ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_MULTIBODY_PARALLEL")
    _MULTIBODY_MODE_CACHE[] = mode
    return mode
end

@inline function _multibody_parallel_mode()::Symbol
    cached = _MULTIBODY_MODE_CACHE[]
    return cached === nothing ? refresh_multibody_parallel_mode!() : cached
end


@inline function _multibody_thread_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD", 4)
end

@inline function _multibody_max_threads()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_MAX_THREADS", 4)
end

@inline function _threadid_capacity()::Int
    # Compatibility shim for legacy tests/callers after migrating scratch indexing
    # from thread ids to stable worker ids.
    return max(Threads.maxthreadid(), _multibody_max_threads())
end

@inline function _multibody_outer_parallel_hint()::Bool
    # Explicit hint for disabling intra-satellite threading under outer parallel execution.
    return ParallelPolicy.outer_parallel_active()
end

@inline function _multibody_use_threads(num_items::Int; heavy_work::Bool=true)::Bool
    return _multibody_thread_decision(num_items; heavy_work=heavy_work).use_threads
end

@inline function _multibody_thread_decision(num_items::Int; heavy_work::Bool=true)
    mode = _multibody_parallel_mode()

    # Cheap forced-false checks before the policy call, because this runs PER
    # SATELLITE PER RHS CALL.
    #
    # Both call sites -- the n-body force and the aerodynamic wrench -- invoke
    # this from inside the per-satellite effector chain, and the answer does not
    # depend on the satellite: it is a function of `num_items` (third bodies, or
    # aerodynamic links), which is the same for every one of them. Reaching
    # thread_policy_decision costs roughly nine ENV reads, each allocating
    # through strip/lowercase, plus two telemetry lock acquisitions. On a
    # 256-satellite constellation with third-body gravity that is ~2300 ENV
    # reads and 512 lock acquisitions per RHS call, to re-derive one answer 256
    # times.
    #
    # The common case never had a decision to make. SPACEAGORA_MULTIBODY_THREAD_THRESHOLD
    # defaults to 4 and a typical third-body set is Sun plus Moon, so num_items
    # is 2 and threading is refused on the threshold -- after all that work.
    #
    # `mode == :on` is excluded from the threshold check because it forces
    # threading regardless of item count, so its answer is not forced here.
    if mode == :off || num_items <= 1 || Threads.nthreads() <= 1
        return (use_threads=false, allotment=1, mode=mode)
    end
    threshold = _multibody_thread_threshold()
    if mode != :on && num_items < max(1, threshold)
        return (use_threads=false, allotment=1, mode=mode)
    end
    outer_active = _multibody_outer_parallel_hint()
    allow_with_outer = _parse_bool_env("SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER", false)
    heavy_only = _parse_bool_env("SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY", true)
    policy = ParallelPolicy.thread_policy_decision(
        num_items;
        mode=mode,
        threshold=threshold,
        heavy_work=heavy_work,
        heavy_only=heavy_only,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:multibody
    )
    allotment = min(policy.allotment, _multibody_max_threads())
    use_threads = policy.use_threads && allotment > 1
    return (use_threads=use_threads, allotment=use_threads ? allotment : 1, mode=mode)
end
