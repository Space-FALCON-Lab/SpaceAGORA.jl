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
            error_message="",
            # Taken from the solve's own return value rather than by calling
            # policy_telemetry_snapshot() out here. _active_policy_context()
            # resolves through task-local storage first, and the engine runs the
            # solve under its own scoped PolicyContext, so every routing decision
            # is recorded into that scope -- a snapshot taken from the harness's
            # task afterwards reads the *global* context, which never saw any of
            # it and reports zeros. run_simulation already captures the snapshot
            # from inside the scope when return_solver_metadata=true.
            policy=get(timed.value.result, :parallel_policy, nothing)
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
        error_message=errmsg,
        policy=nothing
    )
end

function ppc_process_sample_task(case_name::String, cfg::PPCConfig, mode_name::String, sample_idx::Int, sample_seed::Int)
    mode = ppc_mode_specs()[mode_name]
    # No outer_tasks override: on a Distributed worker the outer split really is
    # active (this process is one of cfg.process_workers running concurrently),
    # which is exactly what the flag is supposed to say.
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
                # cases.jl's own eager-GRAM-load check inspects ARGS to decide
                # whether this process will run a GRAM-live case (see
                # PPC_GRAM_LIVE_CASES) -- but a Distributed worker started via
                # addprocs doesn't inherit the launching process's command-line
                # ARGS, so that check silently never fires here, and GRAM-live
                # cases dispatched under outer_process crash with `MethodError:
                # no method matching GRAMAtmosphereModel(; planet_name::String)`
                # the same way the world-age bug did before that check existed.
                # Distributed workers are only ever spawned for outer_process
                # batches, so eagerly loading GRAMSuite unconditionally here
                # (rather than threading "does this case need it" through) is
                # cheap relative to the alternative of getting it wrong.
                # @eval, not invokelatest: `Base.invokelatest(f)` only defers
                # the *call*, but evaluating the bare identifier
                # `ppc_ensure_gramsuite_loaded!` to get `f` in the first place
                # still happens in this closure's original (pre-`include`)
                # world, so it throws the same `UndefVarError` before
                # invokelatest ever runs. @eval re-resolves the whole
                # expression, name lookup included, fresh against the current
                # global bindings -- the same reason ppc_ensure_gramsuite_loaded!
                # itself uses `@eval import GRAMSuite` rather than a plain
                # `import`.
                @eval ppc_ensure_gramsuite_loaded!()
                nothing
            end
        end
    end
    return [w for w in workers() if w != myid()]
end

"""
    ppc_resolve_outer_backend(case, cfg, mode, sample_count) -> PPCModeSpec

Resolve a mode whose `backend` is `"auto"` into a concrete backend by asking the
shipped outer-route selector, and return the mode with that answer substituted in.

Modes with an explicit backend are returned unchanged.

This exists because `ppc_run_sample_batch` dispatches the outer split itself, by
branching on `mode.backend`. R4 and R5 were declared here with
`backend="threads"`, which meant the two profiles whose entire purpose is to
*choose* an outer route were hard-pinned to one, and `select_outer_route!` was
never reached -- the harness answered the question it was supposed to be asking.
That is not a subtle bias: on `independent_1sat_1hr` at 256 samples the router
picks `:process` (8.1x) while the pin forced `:threads` (2.0x), and the resulting
"router regret" was a property of this table rather than of the router.

The shipped profile definitions already say `outer_backend=:auto` for R3/R4/R5
(`src/parallel/routing/profile_definitions.jl`), so this restores agreement with
them for the two adaptive profiles. R3 deliberately keeps its explicit `threads`
pin: it is the *static* hybrid baseline the adaptive profiles are measured
against, so it has to stay fixed.

Single-sample batches resolve to `"threads"` unconditionally. There is no outer
split to route when there is one simulation, so the resolution would be moot, and
keeping the literal string identical to the pre-fix value keeps the recorded env
comparable across the constellation phases.
"""
function ppc_resolve_outer_backend(
    case::PPCCaseSpec, cfg::PPCConfig, mode::PPCModeSpec, sample_count::Int
)::PPCModeSpec
    mode.backend == "auto" || return mode
    sample_count > 1 || return _ppc_mode_with_backend(mode, "threads")
    route = try
        probe = ppc_single_config(case.name, cfg; seed=cfg.worker_seed, mc_index=1)
        features = SpaceAGORA.SimulationCampaigns.campaign_route_features(
            probe; samples=sample_count
        )
        PP = SpaceAGORA.ParallelProfiles
        # select_outer_route! rather than default_outer_route: it is the entry
        # point the shipped campaign runner calls, so it carries the candidate
        # filtering and machine clamping too. The state is fresh per worker
        # process, so no cross-run feedback is folded in -- each measured point
        # sees the router's cold decision, which is what makes the points
        # comparable to each other.
        PP.select_outer_route!(
            PP.OuterRouteState(),
            features;
            machine_class=PP._machine_parallel_class(),
            threads_available=Threads.nthreads() > 1,
            parallel_enabled=true,
        )
    catch err
        @warn "Outer-route resolution failed; falling back to threads" case=case.name mode=mode.name exception=err
        :threads
    end
    backend = route === :process ? "process" : (route === :threads ? "threads" : "none")
    return _ppc_mode_with_backend(mode, backend)
end

function _ppc_mode_with_backend(mode::PPCModeSpec, backend::String)::PPCModeSpec
    return PPCModeSpec(;
        (f => (f === :backend ? backend : getfield(mode, f)) for f in fieldnames(PPCModeSpec))...
    )
end

function ppc_run_sample_batch(case::PPCCaseSpec, cfg::PPCConfig, mode::PPCModeSpec, sample_count::Int)
    sample_indices = collect(1:sample_count)
    sample_seeds = [cfg.worker_seed + i - 1 for i in sample_indices]
    # Adaptive profiles (backend="auto") ask the router here; every other mode
    # keeps the backend its spec declares. Done once, before warm-up, so the
    # warm-up, the env recorded in the row, and the timed dispatch below all
    # agree on one backend.
    mode = ppc_resolve_outer_backend(case, cfg, mode, sample_count)
    uses_process_pool = sample_count > 1 && mode.backend == "process"
    # How many outer units of work this batch actually dispatches concurrently,
    # mirroring the branch conditions below exactly. Single-sample batches --
    # which is every constellation case in the catalog -- run one simulation
    # with nothing beside it, so no outer split is active no matter which mode
    # is nominally under test. See ppc_mode_env_pairs for why that distinction
    # changes the measurement.
    outer_tasks = if uses_process_pool
        sample_count
    elseif sample_count > 1 && mode.backend == "threads" && Threads.nthreads() > 1
        sample_count
    else
        1
    end
    withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=outer_tasks)...) do
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
    # Provision and warm the Distributed pool *before* the clock starts. Both
    # costs are one-off per worker process and neither is part of the throughput
    # this phase claims to measure: addprocs plus the workers' `include` of the
    # study files and `using SpaceAGORA` runs tens of seconds, and each worker
    # then JIT-compiles the whole RHS/solver stack on its first solve. Left
    # inside the timed region they dominate every low-sample-count point and
    # make the throughput-vs-workers curve bend the wrong way (more workers =
    # more startup to pay for).
    if uses_process_pool
        pre_workers = ppc_ensure_process_workers!(cfg.process_workers)
        if cfg.warmup > 0 && !isempty(pre_workers)
            @sync for (offset, w) in enumerate(pre_workers)
                @async try
                    remotecall_wait(
                        ppc_process_sample_task,
                        w,
                        case.name,
                        cfg,
                        mode.name,
                        offset,
                        cfg.worker_seed - offset,
                    )
                catch
                end
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
        withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=outer_tasks)...) do
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
            withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=outer_tasks)...) do
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
        withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=outer_tasks)...) do
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
        execution_scope=execution_scope,
        outer_tasks=outer_tasks,
        policy=ppc_policy_columns(
            isempty(results) ? nothing :
                (results[end] isa NamedTuple && haskey(results[end], :policy) ? results[end].policy : nothing)
        )
    )
end

# What the router actually decided during the timed batch.
#
# Until this existed the CSV recorded `rhs_execution_mode="auto"` -- the mode
# *requested*, never the route selected -- so a finding like review point 8's "it
# selects poorly for the GRAM case" could be measured but not attributed. The
# selected RHS plan lives in per-solve `shared_buffers.rhs_plan_override`, which
# is gone by the time the solve returns; `record_rhs_plan_selection!`
# (src/parallel/policy/policy_telemetry.jl) copies it into the process-level
# policy telemetry so it survives, and this reads it back.
#
# Process-route caveat: the counters are per-process, and on the `outer_process`
# route the solves happen on Distributed workers, so the controller's snapshot is
# empty and these columns come back as the `:none`/0 defaults. That is correct
# rather than missing -- the controller process genuinely made no routing
# decision -- but it means the selected-route columns are only meaningful for the
# in-process (serial/threads) routes.
function ppc_policy_columns(snap)
    empty = (
        rhs_plan_source="none", rhs_plan_mode="none", rhs_plan_allotment=0,
        policy_last_mode="none", policy_last_allotment=0, policy_last_outer_active=false,
        policy_decisions_total=0, policy_adaptive_decisions_total=0,
        policy_threads_enabled_total=0, policy_discarded_by_route_total=0,
    )
    snap === nothing && return empty
    g(k, d) = hasproperty(snap, k) ? getproperty(snap, k) : d
    return (
        rhs_plan_source=string(g(:rhs_plan_source, "none")),
        rhs_plan_mode=string(g(:rhs_plan_mode, "none")),
        rhs_plan_allotment=g(:rhs_plan_allotment, 0),
        policy_last_mode=string(g(:last_mode, "none")),
        policy_last_allotment=g(:last_allotment, 0),
        policy_last_outer_active=g(:last_outer_active, false),
        policy_decisions_total=g(:decisions_total, 0),
        policy_adaptive_decisions_total=g(:adaptive_decisions_total, 0),
        policy_threads_enabled_total=g(:threads_enabled_total, 0),
        policy_discarded_by_route_total=g(:policy_discarded_by_route_total, 0),
    )
end

function ppc_run_worker_performance(cfg::PPCConfig)
    catalog = ppc_case_catalog()
    case = catalog[cfg.worker_case]
    mode = ppc_mode_specs()[cfg.worker_mode]
    rows = NamedTuple[]
    hw = ppc_hardware_snapshot()
    samples = case.montecarlo ? cfg.worker_mc_samples : 1

    # One row per timed repeat, all inside this one process. See
    # PPCConfig.worker_repeats for why the repeats are not separate subprocesses.
    # ppc_run_sample_batch re-runs its own warm-up on each pass, so every repeat
    # is measured the same way; only the first pays cold JIT, which is exactly the
    # cost this loop exists to stop paying N times.
    for repeat_offset in 0:(max(1, cfg.worker_repeats) - 1)
    repeat_index = cfg.worker_repeat + repeat_offset
    # Warm-up runs on the first repeat only. Its job is to get the RHS/solver
    # stack compiled and the caches populated before anything is timed, and that
    # is a property of the process, not of the repeat -- once repeat 1 has run,
    # repeats 2..N are already warm. Repeating it would just re-pay the solve cost
    # (2 warm-ups x 3 repeats = 6 untimed solves per point instead of 2), which on
    # the heavier cases is a larger cost than the startup this loop exists to
    # amortise.
    repeat_cfg = repeat_offset == 0 ? cfg :
        _ppc_with(cfg;
            worker_repeat=repeat_index,
            worker_seed=cfg.worker_seed + repeat_offset,
            warmup=0,
        )

    batch = ppc_run_sample_batch(case, repeat_cfg, mode, samples)
    # Built from the batch's own outer_tasks so the recorded env matches what
    # the timed run actually saw (see ppc_mode_env_pairs).
    env_string = ppc_effective_env_string(mode, repeat_cfg; outer_tasks=batch.outer_tasks)
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
        repeat=repeat_index,
        seed=repeat_cfg.worker_seed,
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
        # What the router actually selected during the timed batch, as opposed to
        # the mode requested above. See ppc_policy_snapshot for the mechanism and
        # its process-route caveat.
        batch.policy...,
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
    end
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

    withenv(ppc_mode_env_pairs(serial, cfg; outer_tasks=1)...) do
        try
            ref_result = ppc_solve_once(args_ref, cfg)
            ref_ok = ppc_solve_success(ref_result.solution)
            ref_reason = string(ref_result.solution.retcode)
        catch err
            ref_reason = string(typeof(err), ": ", sprint(showerror, err))
        end
    end
    withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=1)...) do
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

function ppc_worker_cmd(cfg::PPCConfig; case::String, mode::String, threads::Int, repeat::Int, seed::Int, mc_samples::Int, outfile::String, parity::Bool, repeats::Int=1)
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
        "--worker-repeats=$(repeats)",
        "--worker-seed=$(seed)",
        "--worker-mc-samples=$(mc_samples)",
        # Must be forwarded explicitly: the worker is a fresh subprocess that
        # re-parses its own CLI from scratch, so anything the controller resolved
        # (here, the phase's warm-up count) is invisible to it unless it appears
        # in argv. Without this the worker fell back to zero warm-up and timed
        # its own JIT compilation instead of the simulation.
        "--warmup=$(cfg.warmup)",
        "--solver-mode=$(cfg.solver_mode)",
        "--process-workers=$(cfg.process_workers)",
        "--parity-samples=$(cfg.parity_samples)",
        "--outfile=$(outfile)",
        "--parity=$(parity ? 1 : 0)"
    ]
    return Cmd(_ppc_apply_cpu_pinning(argv, cfg.cpu_pinning, threads))
end

# Resume support: a worker's outfile is considered already done only if it
# parses and has at least one row with success not missing/false — a file
# left behind by a killed or crashed worker (empty, header-only, or a
# `success=false` row) is treated as not done and re-run, same as if it
# were never there. Any read/parse failure is treated the same way (re-run)
# rather than raising, since this only gates a skip decision.
function _ppc_worker_already_done(outfile::String)::Bool
    isfile(outfile) || return false
    try
        df = ppc_read_optional(outfile)
        hasproperty(df, :success) || return false
        successes = collect(skipmissing(df.success))
        return !isempty(successes) && all(successes)
    catch
        return false
    end
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

    # Refuse to time anything on a machine that is already loaded. See
    # ppc_assert_machine_quiet! -- this is the check whose absence let an entire
    # 115-point run be collected against another user's 26-hour job.
    ppc_assert_machine_quiet!()
    println("[parallelization-performance] profile=$(cfg.profile) outdir=$(outdir)")
    println("[parallelization-performance] load=$(round(ppc_load_average(); digits=2)) " *
            "headroom=$(round(100 * ppc_load_headroom(); digits=0))% " *
            "cores=$(_ppc_physical_core_count())")
    println("[parallelization-performance] cases=$(join(cases, ","))")
    println("[parallelization-performance] modes=$(join(cfg.modes, ","))")

    perf_paths = String[]
    for case in cases
        case_meta = ppc_case_catalog()[case]
        is_mc = case_meta.montecarlo
        # A joint_routing rung is a (spacecraft, samples) pair whose product is
        # the fixed total work the grid compares aspect ratios at, so running it
        # at the profile's sample ladder would measure a different total at every
        # rung and the comparison the grid exists for would not hold. Those cases
        # therefore supply their own sample count, unless the caller named one.
        sample_counts = if !is_mc
            [1]
        elseif case_meta.family == "joint_routing" && !cfg.mc_samples_explicit
            [case_meta.default_samples]
        else
            cfg.mc_samples
        end
        # One subprocess per (case, mode, threads, mc) point, running all of that
        # point's repeats inside it -- not one subprocess per repeat. A worker
        # spends ~80 s on Julia startup plus JIT of the RHS/solver stack before it
        # can time anything, so per-repeat subprocesses paid that cost `repeats`
        # times to collect `repeats` samples of the same point; at three repeats
        # that was two thirds of the entire run. Resume granularity is coarser as a
        # result (a point's repeats complete together or not at all), which is the
        # trade this makes. See PPCConfig.worker_repeats.
        for mode in cfg.modes, thread_count in cfg.threads, mc_count in sample_counts
            if mode == "serial" && thread_count != minimum(cfg.threads)
                continue
            end
            outfile = joinpath(scratch, "perf_$(case)_$(mode)_t$(thread_count)_mc$(mc_count).csv")
            if _ppc_worker_already_done(outfile)
                println("[skip] $(case) mode=$(mode) threads=$(thread_count) mc=$(mc_count) (already completed — resume)")
            else
                cmd = ppc_worker_cmd(
                    cfg;
                    case=case,
                    mode=mode,
                    threads=thread_count,
                    repeat=1,
                    repeats=cfg.repeats,
                    seed=cfg.seed + 1,
                    mc_samples=mc_count,
                    outfile=outfile,
                    parity=false
                )
                println("[run] $(case) mode=$(mode) threads=$(thread_count) repeats=$(cfg.repeats) mc=$(mc_count)")
                run(cmd)
                on_run_complete === nothing || on_run_complete()
            end
            push!(perf_paths, outfile)
        end
    end

    parity_paths = String[]
    parity_modes = [m for m in cfg.modes if m != "serial"]
    for case in parity_cases, mode in parity_modes
        thread_count = maximum(cfg.threads)
        outfile = joinpath(scratch, "parity_$(case)_$(mode)_t$(thread_count).csv")
        if _ppc_worker_already_done(outfile)
            println("[skip] parity $(case) mode=$(mode) threads=$(thread_count) (already completed — resume)")
        else
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
        end
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
