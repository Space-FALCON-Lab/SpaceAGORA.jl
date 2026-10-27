module BenchmarkRouting
export partition_sweep_routes

"""Partition selected cases without executing work or mutating the environment.

Indices preserve the caller's selected order. Unsafe threaded cases fall back to
serial; process requests remain process requests. Execution and reference runs
belong to the individual experiment.
"""
function partition_sweep_routes(selected::AbstractVector; choose_route, thread_safe)
    chosen_routes = fill(:none, length(selected))
    process_tasks = Tuple{Int, eltype(selected)}[]
    thread_indices = Int[]
    serial_indices = Int[]

    for (idx, base_case) in enumerate(selected)
        route = choose_route(base_case)
        chosen_routes[idx] = route
        if route == :process
            push!(process_tasks, (idx, base_case))
        elseif route == :threads
            if thread_safe(base_case)
                push!(thread_indices, idx)
            else
                push!(serial_indices, idx)
                chosen_routes[idx] = :none
            end
        else
            push!(serial_indices, idx)
        end
    end

    return (; chosen_routes, process_tasks, thread_indices, serial_indices)
end
end
