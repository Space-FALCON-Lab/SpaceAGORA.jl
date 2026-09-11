function ppc_solve_once(args, cfg::PPCConfig)
    # The worker applies one mode environment around the whole sample batch.
    # Calling the config-taking API here would reapply those values by mutating
    # process-global ENV for every solve, which is unsafe when samples run on
    # multiple threads.
    return SimulationEngine.run_simulation(
        args;
        isolate_state=false,
        return_solution=true,
        return_solver_metadata=true
    )
end

@inline function ppc_solve_success(sol)::Bool
    return string(sol.retcode) == "Success"
end

function ppc_terminal_metrics(sol)
    if isempty(sol.u)
        return (terminal_time_s=missing, pos_norm_m=missing, vel_norm_mps=missing, mass_kg=missing)
    end
    sc = sol.u[end].sc[1]
    pos = _ppc_vec3(sc, :pos)
    vel = _ppc_vec3(sc, :vel)
    mass = _ppc_mass(sc)
    return (
        terminal_time_s=isempty(sol.t) ? missing : Float64(sol.t[end]),
        pos_norm_m=pos === nothing ? missing : norm(pos),
        vel_norm_mps=vel === nothing ? missing : norm(vel),
        mass_kg=mass
    )
end

function ppc_run_sample_once(case_name::String, cfg::PPCConfig, sample_idx::Int, sample_seed::Int)
    args = ppc_single_config(case_name, cfg; seed=sample_seed, mc_index=sample_idx)
    timed = @timed begin
        try
            result = ppc_solve_once(args, cfg)
            (ok=ppc_solve_success(result.solution), result=result, err=nothing)
        catch err
            (ok=false, result=nothing, err=err)
        end
    end
    if timed.value.ok && timed.value.result !== nothing
        sol = timed.value.result.solution
        return (
            success=true,
            retcode=string(sol.retcode),
            wall_time_s=Float64(timed.time),
            terminal=ppc_terminal_metrics(sol),
            error_type="",
            error_message=""
        )
    end
    retcode = timed.value.result === nothing ? string(typeof(timed.value.err)) : string(timed.value.result.solution.retcode)
    errmsg = timed.value.err === nothing ? "" : sprint(showerror, timed.value.err)
    return (
        success=false,
        retcode=retcode,
        wall_time_s=Float64(timed.time),
        terminal=(terminal_time_s=missing, pos_norm_m=missing, vel_norm_mps=missing, mass_kg=missing),
        error_type=retcode,
        error_message=errmsg
    )
end

function ppc_process_sample_task(case_name::String, cfg::PPCConfig, mode_name::String, sample_idx::Int, sample_seed::Int)
    mode = ppc_mode_specs()[mode_name]
    withenv(ppc_mode_env_pairs(mode, cfg)...) do
        return ppc_run_sample_once(case_name, cfg, sample_idx, sample_seed)
    end
end

function ppc_ensure_process_workers!(n::Int)::Vector{Int}
    desired = max(1, n)
    current_external = max(0, nprocs() - 1)
    if current_external < desired
        add_count = desired - current_external
        new_workers = addprocs(
            add_count;
            exeflags=Cmd(["--threads=1", "--startup-file=no", "--project=$(PPC_REPO_ROOT)"])
        )
        @sync for w in new_workers
            @async remotecall_wait(w, PPC_REPO_ROOT) do repo_root
                study_dir = joinpath(repo_root, "benchmarks", "studies", "parallelization_performance")
                include(joinpath(study_dir, "cli.jl"))
                include(joinpath(study_dir, "modes.jl"))
                include(joinpath(study_dir, "cases.jl"))
                include(joinpath(study_dir, "trajectory_parity.jl"))
                include(joinpath(study_dir, "execution.jl"))
                nothing
            end
        end
    end
    return [w for w in workers() if w != myid()]
end

function ppc_run_sample_batch(case::PPCCaseSpec, cfg::PPCConfig, mode::PPCModeSpec, sample_count::Int)
    sample_indices = collect(1:sample_count)
    sample_seeds = [cfg.worker_seed + i - 1 for i in sample_indices]
    withenv(ppc_mode_env_pairs(mode, cfg)...) do
        for warmup_idx in 1:cfg.warmup
            warmup_args = ppc_single_config(
                case.name,
                cfg;
                seed=cfg.worker_seed - warmup_idx,
                mc_index=warmup_idx
            )
            try
                ppc_solve_once(warmup_args, cfg)
            catch
            end
        end
    end
    GC.gc()
    batch_started = time()
    results = Vector{Any}(undef, sample_count)
    actual_backend = "serial"
    execution_scope = sample_count == 1 ? "single_simulation" : "serial_sample_batch"

    if sample_count > 1 && mode.backend == "threads" && Threads.nthreads() > 1
        actual_backend = "threads"
        execution_scope = "outer_thread_sample_batch"
        withenv(ppc_mode_env_pairs(mode, cfg)...) do
            Threads.@threads for i in eachindex(sample_indices)
                results[i] = ppc_run_sample_once(case.name, cfg, sample_indices[i], sample_seeds[i])
            end
        end
    elseif sample_count > 1 && mode.backend == "process"
        actual_backend = "process"
        execution_scope = "outer_process_sample_batch"
        worker_ids = ppc_ensure_process_workers!(cfg.process_workers)
        if isempty(worker_ids)
            actual_backend = "serial"
            execution_scope = "serial_sample_batch_process_unavailable"
            withenv(ppc_mode_env_pairs(mode, cfg)...) do
                for i in eachindex(sample_indices)
                    results[i] = ppc_run_sample_once(case.name, cfg, sample_indices[i], sample_seeds[i])
                end
            end
        else
            tasks = [(case.name, cfg, mode.name, sample_indices[i], sample_seeds[i]) for i in eachindex(sample_indices)]
            pool = CachingPool(worker_ids)
            results = pmap(pool, tasks) do task
                case_name, local_cfg, mode_name, sample_idx, sample_seed = task
                ppc_process_sample_task(case_name, local_cfg, mode_name, sample_idx, sample_seed)
            end
        end
    else
        withenv(ppc_mode_env_pairs(mode, cfg)...) do
            for i in eachindex(sample_indices)
                results[i] = ppc_run_sample_once(case.name, cfg, sample_indices[i], sample_seeds[i])
            end
        end
    end

    batch_wall = Float64(time() - batch_started)
    return (
        results=results,
        batch_wall_time_s=batch_wall,
        actual_backend=actual_backend,
        execution_scope=execution_scope
    )
end

function ppc_run_worker_performance(cfg::PPCConfig)
    catalog = ppc_case_catalog()
    case = catalog[cfg.worker_case]
    mode = ppc_mode_specs()[cfg.worker_mode]
    rows = NamedTuple[]
    hw = ppc_hardware_snapshot()
    env_string = ppc_effective_env_string(mode, cfg)
    samples = case.montecarlo ? cfg.worker_mc_samples : 1

    batch = ppc_run_sample_batch(case, cfg, mode, samples)
    sample_results = batch.results
    total_success = all(r -> r.success, sample_results)
    sample_wall_sum = sum(r -> Float64(r.wall_time_s), sample_results)
    final_result = sample_results[end]
    final_retcode = total_success ? string(final_result.retcode) : join(unique(string(r.retcode) for r in sample_results if !r.success), "|")
    final_terminal = final_result.terminal
    throughput = samples / max(batch.batch_wall_time_s, eps(Float64))

    push!(rows, (
        timestamp_utc=hw.timestamp_utc,
        git_commit=hw.git_commit,
        machine=hw.machine,
        julia_version=hw.julia_version,
        cpu_threads=hw.cpu_threads,
        julia_threads=hw.julia_threads,
        case=cfg.worker_case,
        family=case.family,
        mode=cfg.worker_mode,
        parallel_profile=mode.profile,
        thread_count=cfg.worker_threads,
        process_workers=cfg.process_workers,
        repeat=cfg.worker_repeat,
        seed=cfg.worker_seed,
        mc_samples=samples,
        solver_mode=cfg.solver_mode,
        success=total_success,
        retcode=final_retcode,
        wall_time_s=batch.batch_wall_time_s,
        sample_wall_time_sum_s=sample_wall_sum,
        mean_sample_wall_time_s=sample_wall_sum / max(1, samples),
        execution_scope=batch.execution_scope,
        outer_backend_actual=batch.actual_backend,
        outer_tasks=samples,
        throughput_samples_per_s=throughput,
        rhs_execution_mode="auto",
        rhs_batch_parallel=mode.rhs_batch,
        density_callback_parallel=mode.density,
        control_callback_parallel=mode.control,
        thermal_callback_parallel=mode.thermal,
        effector_parallel=mode.effector,
        inner_allowed_with_outer=mode.allow_inner_with_outer,
        terminal_time_s=final_terminal.terminal_time_s,
        final_primary_pos_norm_m=final_terminal.pos_norm_m,
        final_primary_vel_norm_mps=final_terminal.vel_norm_mps,
        final_primary_mass_kg=final_terminal.mass_kg,
        effective_env=env_string
    ))
    ppc_write_rows(cfg.worker_outfile, rows)
    return nothing
end

function ppc_run_worker_parity(cfg::PPCConfig)
    mode = ppc_mode_specs()[cfg.worker_mode]
    serial = ppc_mode_specs()["serial"]
    case = ppc_case_catalog()[cfg.worker_case]
    args_ref = ppc_single_config(cfg.worker_case, cfg; seed=cfg.worker_seed, mc_index=1)
    args_cmp = deepcopy(args_ref)
    ref_result = nothing
    cmp_result = nothing
    ref_ok = false
    cmp_ok = false
    ref_reason = ""
    cmp_reason = ""

    withenv(ppc_mode_env_pairs(serial, cfg)...) do
        try
            ref_result = ppc_solve_once(args_ref, cfg)
            ref_ok = ppc_solve_success(ref_result.solution)
            ref_reason = string(ref_result.solution.retcode)
        catch err
            ref_reason = string(typeof(err), ": ", sprint(showerror, err))
        end
    end
    withenv(ppc_mode_env_pairs(mode, cfg)...) do
        try
            cmp_result = ppc_solve_once(args_cmp, cfg)
            cmp_ok = ppc_solve_success(cmp_result.solution)
            cmp_reason = string(cmp_result.solution.retcode)
        catch err
            cmp_reason = string(typeof(err), ": ", sprint(showerror, err))
        end
    end

    rows = NamedTuple[]
    if ref_ok && cmp_ok
        metrics = ppc_compare_trajectories(ref_result.solution, cmp_result.solution, args_ref; sample_count=cfg.parity_samples)
        push!(rows, merge((
            case=cfg.worker_case,
            family=case.family,
            mode=cfg.worker_mode,
            thread_count=cfg.worker_threads,
            seed=cfg.worker_seed,
            solver_mode=cfg.solver_mode,
            reference_retcode=ref_reason,
            candidate_retcode=cmp_reason,
        ), metrics))
    else
        push!(rows, (
            case=cfg.worker_case,
            family=case.family,
            mode=cfg.worker_mode,
            thread_count=cfg.worker_threads,
            seed=cfg.worker_seed,
            solver_mode=cfg.solver_mode,
            reference_retcode=ref_reason,
            candidate_retcode=cmp_reason,
            pass=false,
            samples=0,
            pos_rel_rms=missing,
            pos_rel_p90=missing,
            pos_rel_max=missing,
            vel_rel_rms=missing,
            vel_rel_p90=missing,
            vel_rel_max=missing,
            q_angle_max_rad=missing,
            omega_rel_max=missing,
            mass_rel_max=missing,
            event_count_equal=false,
            event_time_abs_max_s=missing,
            ref_periapsis_count=missing,
            cmp_periapsis_count=missing,
            ref_interface_count=missing,
            cmp_interface_count=missing
        ))
    end
    ppc_write_rows(cfg.worker_outfile, rows)
    return nothing
end

# Prefixes a worker's argv with `taskset -c <cpus>` when a CPU pool was
# reserved via --cpu-list/SPACEAGORA_PPC_CPU_LIST, pinning that worker to the
# first `threads` cores of the pool. Linux-only (taskset); pinning is skipped
# with a warning on other platforms or if taskset isn't on PATH.
function _ppc_apply_cpu_pinning(argv::Vector{String}, cpu_pinning::Vector{Int}, threads::Int)::Vector{String}
    isempty(cpu_pinning) && return argv
    if !Sys.islinux()
        @warn "CPU pinning requested but taskset is only available on Linux; running unpinned."
        return argv
    end
    if Sys.which("taskset") === nothing
        @warn "CPU pinning requested but taskset was not found on PATH; running unpinned."
        return argv
    end
    threads <= length(cpu_pinning) || throw(ArgumentError(
        "CPU pinning pool has $(length(cpu_pinning)) core(s) but this run needs $(threads) thread(s); " *
        "pass a larger --cpu-list, or drop it to disable pinning."
    ))
    cpu_list = join(cpu_pinning[1:threads], ",")
    return vcat(["taskset", "-c", cpu_list], argv)
end

function ppc_worker_cmd(cfg::PPCConfig; case::String, mode::String, threads::Int, repeat::Int, seed::Int, mc_samples::Int, outfile::String, parity::Bool)
    julia_bin = Base.julia_cmd().exec[1]
    argv = String[
        julia_bin,
        "--threads=$(threads)",
        "--project=$(PPC_REPO_ROOT)",
        PPC_LAUNCHER,
        cfg.profile,
        "--worker",
        "--case=$(case)",
        "--mode=$(mode)",
        "--thread-count=$(threads)",
        "--repeat=$(repeat)",
        "--worker-seed=$(seed)",
        "--worker-mc-samples=$(mc_samples)",
        "--solver-mode=$(cfg.solver_mode)",
        "--process-workers=$(cfg.process_workers)",
        "--parity-samples=$(cfg.parity_samples)",
        "--outfile=$(outfile)",
        "--parity=$(parity ? 1 : 0)"
    ]
    return Cmd(_ppc_apply_cpu_pinning(argv, cfg.cpu_pinning, threads))
end

function ppc_run_controller(cfg::PPCConfig; on_run_complete::Union{Nothing, Function}=nothing)
    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    outdir = cfg.outdir == PPC_DEFAULT_OUTDIR ? joinpath(PPC_DEFAULT_OUTDIR, stamp) : cfg.outdir
    mkpath(outdir)
    scratch = joinpath(outdir, "worker_rows")
    mkpath(scratch)
    cases = ppc_resolve_cases(cfg.cases)
    parity_cases = ppc_resolve_cases(cfg.parity_cases)
    modes = ppc_mode_specs()
    unknown_modes = [m for m in cfg.modes if !haskey(modes, m)]
    isempty(unknown_modes) || throw(ArgumentError("Unknown mode(s): $(join(unknown_modes, ", "))"))

    println("[parallelization-performance] profile=$(cfg.profile) outdir=$(outdir)")
    println("[parallelization-performance] cases=$(join(cases, ","))")
    println("[parallelization-performance] modes=$(join(cfg.modes, ","))")

    perf_paths = String[]
    for case in cases
        is_mc = ppc_case_catalog()[case].montecarlo
        sample_counts = is_mc ? cfg.mc_samples : [1]
        for mode in cfg.modes, thread_count in cfg.threads, mc_count in sample_counts, repeat in 1:cfg.repeats
            if mode == "serial" && thread_count != minimum(cfg.threads)
                continue
            end
            outfile = joinpath(scratch, "perf_$(case)_$(mode)_t$(thread_count)_mc$(mc_count)_r$(repeat).csv")
            cmd = ppc_worker_cmd(
                cfg;
                case=case,
                mode=mode,
                threads=thread_count,
                repeat=repeat,
                seed=cfg.seed + repeat,
                mc_samples=mc_count,
                outfile=outfile,
                parity=false
            )
            println("[run] $(case) mode=$(mode) threads=$(thread_count) repeat=$(repeat) mc=$(mc_count)")
            run(cmd)
            push!(perf_paths, outfile)
            on_run_complete === nothing || on_run_complete()
        end
    end

    parity_paths = String[]
    parity_modes = [m for m in cfg.modes if m != "serial"]
    for case in parity_cases, mode in parity_modes
        thread_count = maximum(cfg.threads)
        outfile = joinpath(scratch, "parity_$(case)_$(mode)_t$(thread_count).csv")
        cmd = ppc_worker_cmd(
            cfg;
            case=case,
            mode=mode,
            threads=thread_count,
            repeat=1,
            seed=cfg.seed,
            mc_samples=1,
            outfile=outfile,
            parity=true
        )
        println("[parity] $(case) mode=$(mode) threads=$(thread_count)")
        run(cmd)
        push!(parity_paths, outfile)
    end

    raw = DataFrame()
    for path in perf_paths
        isfile(path) || continue
        raw = vcat(raw, ppc_read_optional(path); cols=:union)
    end
    parity = DataFrame()
    for path in parity_paths
        isfile(path) || continue
        parity = vcat(parity, ppc_read_optional(path); cols=:union)
    end
    summary = ppc_summarize(raw, parity)
    raw_path = joinpath(outdir, "parallelization_performance_raw_$(cfg.profile)_$(stamp).csv")
    summary_path = joinpath(outdir, "parallelization_performance_summary_$(cfg.profile)_$(stamp).csv")
    parity_path = joinpath(outdir, "parallelization_trajectory_parity_$(cfg.profile)_$(stamp).csv")
    hardware_path = joinpath(outdir, "parallelization_hardware_$(cfg.profile)_$(stamp).csv")
    report_path = joinpath(outdir, "parallelization_performance_report_$(cfg.profile)_$(stamp).md")
    CSV.write(raw_path, raw)
    CSV.write(summary_path, summary)
    CSV.write(parity_path, parity)
    CSV.write(hardware_path, DataFrame([ppc_hardware_snapshot()]))
    plot_paths = ppc_write_plots(outdir, summary, parity)
    ppc_write_report(report_path, cfg, raw, summary, parity)
    println("[done] raw=$(raw_path)")
    println("[done] summary=$(summary_path)")
    println("[done] parity=$(parity_path)")
    println("[done] report=$(report_path)")
    println("[done] plots=$(length(plot_paths))")
    return nothing
end

function main_parallelization_performance()
    cfg = parse_parallelization_performance_cli()
    if cfg.worker
        if isempty(cfg.worker_outfile) || isempty(cfg.worker_case) || isempty(cfg.worker_mode)
            throw(ArgumentError("--worker requires --case, --mode, and --outfile."))
        end
        if cfg.worker_parity
            ppc_run_worker_parity(cfg)
        else
            ppc_run_worker_performance(cfg)
        end
    else
        ppc_run_controller(cfg)
    end
    return nothing
end
