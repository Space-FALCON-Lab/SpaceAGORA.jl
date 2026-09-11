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

function _persistent_worker_pool_for(source::Symbol)::_PersistentForeachPool
    key = _persistent_pool_key(source)
    lock(_persistent_foreach_lock) do
        return get!(_persistent_foreach_worker_pools, key) do
            _create_persistent_foreach_worker_pool(_default_thread_pool_size())
        end
    end
end

function threaded_foreach_worker_persistent(
    source::Symbol,
    num_items::Int,
    allotment::Int,
    f::F
) where {F <: Function}
    num_items <= 0 && return nothing
    workers = _thread_worker_count(num_items, allotment)
    if workers <= 1 || !callback_persistent_workers_enabled()
        return threaded_foreach_worker(num_items, allotment, f)
    end
    pool = _persistent_worker_pool_for(source)
    return _threaded_foreach_persistent!(pool, num_items, workers, f)
end

function threaded_foreach_worker_persistent(
    f::F,
    source::Symbol,
    num_items::Int,
    allotment::Int
) where {F <: Function}
    return threaded_foreach_worker_persistent(source, num_items, allotment, f)
end

# Spin-barrier variant: workers spin-poll an atomic generation counter instead of
# sleeping on a Channel. Dispatch overhead is ~10-50 ns vs ~1-5 µs for channels,
# allowing the harmonics SIMD batch to scale to 32-128+ threads.
# Opt-in via SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER=1.
@inline harmonics_batch_spin_barrier_enabled()::Bool =
    parse_bool_env("SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER", false)

function threaded_foreach_worker_spin(
    source::Symbol,
    num_items::Int,
    allotment::Int,
    f::F,
) where {F <: Function}
    num_items <= 0 && return nothing
    workers = _thread_worker_count(num_items, allotment)
    if workers <= 1
        @inbounds for idx in 1:num_items
            f(1, idx)
        end
        return nothing
    end
    pool = _spin_barrier_pool_for(source)
    return _spin_barrier_dispatch!(pool, num_items, workers, f)
end

function threaded_foreach_worker_spin(
    f::F,
    source::Symbol,
    num_items::Int,
    allotment::Int,
) where {F <: Function}
    return threaded_foreach_worker_spin(source, num_items, allotment, f)
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
