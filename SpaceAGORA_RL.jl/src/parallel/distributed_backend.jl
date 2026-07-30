function setup_distributed_workers(n_workers::Integer)
    n = max(1, Int(n_workers))
    if n == 1
        return workers()
    end
    needed = n - length(workers())
    needed > 0 && addprocs(needed)
    return workers()[1:n]
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
        for process_id in process_ids
            remotecall_wait(Base.eval, process_id, Main, :(using SpaceAGORA_RL))
        end
    catch
        rmprocs(process_ids)
        rethrow()
    end
    return process_ids
end
