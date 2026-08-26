function setup_distributed_workers(n_workers::Integer)
    n = max(1, Int(n_workers))
    if n == 1
        return workers()
    end
    needed = n - length(workers())
    needed > 0 && addprocs(needed)
    return workers()[1:n]
end

function _remotecall_wait_all(f, process_ids::AbstractVector{<:Integer}, args...)
    futures = map(process_ids) do process_id
        remotecall(f, process_id, args...)
    end
    foreach(fetch, futures)
    return nothing
end

function setup_isolated_process_workers(n_workers::Integer)
    n = max(1, Int(n_workers))
    project_file = Base.active_project()
    project_dir = project_file === nothing ? package_root() : dirname(project_file)
    process_ids = addprocs(
        n;
        exeflags=Cmd(["--project=$(project_dir)", "--threads=1"]),
    )
    try
        _remotecall_wait_all(Base.eval, process_ids, Main, :(using SpaceAGORA_RL))
    catch
        rmprocs(process_ids)
        rethrow()
    end
    return process_ids
end
