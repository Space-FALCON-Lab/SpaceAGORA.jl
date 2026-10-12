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

function _persistent_foreach_worker_loop_w(
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
                        f(worker_id, idx)
                    end
                end
            else
                @inbounds for idx in worker_id:active_workers:num_items
                    f(worker_id, idx)
                end
            end
        catch err
            captured = Base.CapturedException(err, catch_backtrace())
        end
        put!(done_channel, captured)
    end
end

function _create_persistent_foreach_worker_pool(workers::Int)::_PersistentForeachPool
    workers = max(2, workers)
    request_channels = [Channel{Any}(1) for _ in 1:workers]
    done_channel = Channel{Any}(workers)
    pool = _PersistentForeachPool(
        workers=workers,
        request_channels=request_channels,
        done_channel=done_channel
    )
    @inbounds for worker_id in 1:workers
        Threads.@spawn _persistent_foreach_worker_loop_w(
            worker_id,
            request_channels[worker_id],
            done_channel
        )
    end
    return pool
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

# ---------------------------------------------------------------------------
# Spin-barrier pool: worker loop, creation, dispatch, teardown.
# ---------------------------------------------------------------------------

function _spin_barrier_worker_loop_w(worker_id::Int, pool::_SpinBarrierPool)::Nothing
    my_gen = 0
    while true
        # Spin-poll until the coordinator bumps our generation counter.
        while pool.worker_gen[worker_id][] == my_gen
            pool.stop[] && return nothing
            GC.safepoint()
        end
        my_gen += 1
        pool.stop[] && return nothing

        request = pool.request[]
        captured = nothing
        try
            num_items      = request.num_items
            active_workers = request.active_workers
            scheduler      = request.scheduler
            chunk          = request.chunk
            f              = request.f
            if scheduler == :dynamic
                next_index = request.next_index
                @inbounds while true
                    start_idx = Threads.atomic_add!(next_index, chunk)
                    start_idx > num_items && break
                    stop_idx = min(num_items, start_idx + chunk - 1)
                    for idx in start_idx:stop_idx
                        f(worker_id, idx)
                    end
                end
            else
                @inbounds for idx in worker_id:active_workers:num_items
                    f(worker_id, idx)
                end
            end
        catch err
            captured = Base.CapturedException(err, catch_backtrace())
        end
        # Record the round's outcome in this worker's private slot (clearing any
        # stale error), then signal done. The atomic_add! release-publishes the
        # slot write; the coordinator reads slots only after the barrier.
        pool.errors[worker_id] = captured
        Threads.atomic_add!(pool.done_count, 1)
    end
end

function _create_spin_barrier_pool(workers::Int)::_SpinBarrierPool
    # Spin workers busy-wait without yielding, so the pool must leave one Julia
    # thread free for the coordinator; more spinners than that would deadlock.
    workers = max(1, min(workers, Threads.nthreads() - 1))
    pool = _SpinBarrierPool(workers)
    @inbounds for worker_id in 1:workers
        Threads.@spawn _spin_barrier_worker_loop_w(worker_id, pool)
    end
    return pool
end

function _spin_barrier_dispatch!(
    pool::_SpinBarrierPool,
    num_items::Int,
    workers::Int,
    f::F,
) where {F <: Function}
    # The coordinator executes the last worker slot itself; only workers-1 pool
    # tasks are signalled. This ensures the coordinator's Julia thread is never
    # blocked waiting for a pool task that has no thread to run on.
    pool_workers = workers - 1
    scheduler  = inner_scheduler_mode()
    chunk      = inner_dynamic_chunk_size()
    next_index = Threads.Atomic{Int}(1)
    lock(pool.run_lock) do
        pool.request[] = (
            num_items      = num_items,
            active_workers = workers,
            scheduler      = scheduler,
            chunk          = chunk,
            next_index     = next_index,
            f              = f,
        )
        # Signal pool workers to start.
        @inbounds for w in 1:pool_workers
            Threads.atomic_add!(pool.worker_gen[w], 1)
        end
        # Coordinator executes the last worker slot (worker_id = workers).
        coordinator_error = nothing
        try
            if scheduler == :dynamic
                @inbounds while true
                    start_idx = Threads.atomic_add!(next_index, chunk)
                    start_idx > num_items && break
                    stop_idx = min(num_items, start_idx + chunk - 1)
                    for idx in start_idx:stop_idx
                        f(workers, idx)
                    end
                end
            else
                @inbounds for idx in workers:workers:num_items
                    f(workers, idx)
                end
            end
        catch err
            coordinator_error = err
        end
        # Spin-wait for pool workers to finish.
        while pool.done_count[] < pool_workers
            GC.safepoint()
        end
        Threads.atomic_sub!(pool.done_count, pool_workers)
        # Propagate the first pool-worker error (by worker index), then any
        # coordinator error. Each slot is written by exactly one worker.
        @inbounds for w in 1:pool_workers
            captured = pool.errors[w]
            if captured !== nothing
                throw(captured.ex)
            end
        end
        if coordinator_error !== nothing
            throw(coordinator_error)
        end
    end
    return nothing
end

function _shutdown_spin_barrier_pool!(pool::_SpinBarrierPool)::Nothing
    pool.stop[] = true
    # Bump all worker generations so they wake from their spin-poll and see stop.
    @inbounds for gen in pool.worker_gen
        Threads.atomic_add!(gen, 1)
    end
    return nothing
end

function _spin_barrier_pool_for(source::Symbol)::_SpinBarrierPool
    key = _persistent_pool_key(source)
    lock(_spin_barrier_lock) do
        return get!(_spin_barrier_pools, key) do
            # Pool needs nthreads-1 workers: the coordinator occupies 1 thread itself,
            # so only nthreads-1 threads are available for spinning pool tasks.
            _create_spin_barrier_pool(max(1, Threads.nthreads() - 1))
        end
    end
end

function _destroy_persistent_foreach_scope!(scope_id::UInt)::Nothing
    channel_pools = _PersistentForeachPool[]
    spin_pools    = _SpinBarrierPool[]
    lock(_persistent_foreach_lock) do
        for dict in (_persistent_foreach_pools, _persistent_foreach_worker_pools)
            stale_keys = Tuple{UInt, Symbol}[]
            for (key, pool) in dict
                if key[1] == scope_id
                    push!(stale_keys, key)
                    push!(channel_pools, pool)
                end
            end
            @inbounds for key in stale_keys
                delete!(dict, key)
            end
        end
    end
    lock(_spin_barrier_lock) do
        stale_keys = Tuple{UInt, Symbol}[]
        for (key, pool) in _spin_barrier_pools
            if key[1] == scope_id
                push!(stale_keys, key)
                push!(spin_pools, pool)
            end
        end
        @inbounds for key in stale_keys
            delete!(_spin_barrier_pools, key)
        end
    end
    @inbounds for pool in channel_pools
        _shutdown_persistent_foreach_pool!(pool)
    end
    @inbounds for pool in spin_pools
        _shutdown_spin_barrier_pool!(pool)
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
