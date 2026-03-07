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
