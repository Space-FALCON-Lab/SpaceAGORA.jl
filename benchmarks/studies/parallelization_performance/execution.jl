function ppc_solve_once(args, cfg::PPCConfig)
    # The worker applies one mode environment around the whole sample batch.
    # Calling the config-taking API here would reapply those values by mutating
    # process-global ENV for every solve, which is unsafe when samples run on
    # multiple threads.
    result = SimulationEngine.run_simulation(
        args;
        isolate_state=false,
        return_solution=true,
        return_solver_metadata=true
    )
    ppc_startup_trace("first_solve")
    return result
end

# SPACEAGORA_PPC_DUMP_STATE_DIR=<dir>: write each successful sample's step times
# and final state, as raw Float64, to <dir>. For byte-for-byte comparison of two
# runs of the same point (the precompile workload's validation compares a run
# with the workload against one without). Off by default.
function _ppc_dump_state(case_name::String, sample_idx::Int, sample_seed::Int, sol)
    dir = get(ENV, "SPACEAGORA_PPC_DUMP_STATE_DIR", "")
    isempty(dir) && return nothing
    mkpath(dir)
    open(joinpath(dir, "state_$(case_name)_i$(sample_idx)_s$(sample_seed).bin"), "w") do io
        write(io, Float64.(sol.t))
        isempty(sol.u) || write(io, Float64[x for x in sol.u[end]])
    end
    return nothing
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
        _ppc_dump_state(case_name, sample_idx, sample_seed, sol)
        return (
            success=true,
            retcode=string(sol.retcode),
            wall_time_s=Float64(timed.time),
            gc_time_s=Float64(timed.gctime),
            alloc_bytes=Int(timed.bytes),
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
        gc_time_s=Float64(timed.gctime),
        alloc_bytes=Int(timed.bytes),
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

# This is a second, independent worker-spawn path from
# SpaceAGORA.ParallelProcess.ensure_process_workers! (used via
# adopt_process_workers! below the outer_process route also touches), with
# its own hardcoded exeflags -- so the pool-worker heap-size-hint default
# has to be applied here too, not just in src/parallel/process/worker_pool.jl.
# It reuses that module's exact formula/constants and its env override
# (SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT), rather than a third knob, since
# these are the same kind of worker under the same per-process memory-share
# reasoning; see docs/architecture/heap_contention.md.
@inline function _ppc_pool_worker_exeflags(desired::Int)::Cmd
    exeflags = String["--threads=1", "--startup-file=no", "--project=$(PPC_REPO_ROOT)"]
    hint = SpaceAGORA.ParallelProcess._pool_worker_heap_size_hint(desired)
    hint === nothing || push!(exeflags, "--heap-size-hint=$(hint)")
    return Cmd(exeflags)
end

function ppc_ensure_process_workers!(n::Int)::Vector{Int}
    desired = max(1, n)
    current_external = max(0, nprocs() - 1)
    if current_external < desired
        add_count = desired - current_external
        new_workers = addprocs(add_count; exeflags=_ppc_pool_worker_exeflags(desired))
        # GRAMSuite goes onto the pool exactly when this process has it, i.e.
        # when this process's case builds a live GRAM atmosphere: cases.jl loads
        # it at include time from the --case on the command line, which a
        # Distributed worker does not see (addprocs does not pass ARGS on), so
        # the answer has to be carried over from here.
        #
        # This used to load it on every pool worker unconditionally, on the
        # reasoning that it was cheap. It is not cheap once the worker has
        # precompiled solver code to lose: the GRAM extension's methods
        # invalidate SpaceAGORA's compiled solve path, and with the paper
        # harness's precompile workload loaded a pool worker's first
        # non-GRAM sample took 19.9 s with GRAMSuite loaded against 1.9 s
        # without (independent_1sat_1hr, one worker, this workstation). Only
        # startup changes; a worker whose case never builds a GRAM model never
        # calls into GRAMSuite.
        load_gram = Base.get_extension(SpaceAGORA, :SpaceAGORAGRAMSuiteExt) !== nothing
        @sync for w in new_workers
            @async remotecall_wait(w, PPC_REPO_ROOT, load_gram) do repo_root, load_gram
                study_dir = joinpath(repo_root, "benchmarks", "studies", "parallelization_performance")
                include(joinpath(study_dir, "cli.jl"))
                include(joinpath(study_dir, "modes.jl"))
                include(joinpath(study_dir, "cases.jl"))
                include(joinpath(study_dir, "trajectory_parity.jl"))
                include(joinpath(study_dir, "execution.jl"))
                # A GRAM-live case dispatched under outer_process needs GRAMSuite
                # on the worker, or it crashes with `MethodError: no method
                # matching GRAMAtmosphereModel(; planet_name::String)` (the
                # extension never attached). @eval, not invokelatest:
                # `Base.invokelatest(f)` only defers the *call*, but evaluating
                # the bare identifier `ppc_ensure_gramsuite_loaded!` to get `f`
                # in the first place still happens in this closure's original
                # (pre-`include`) world, so it throws the same `UndefVarError`
                # before invokelatest ever runs. @eval re-resolves the whole
                # expression, name lookup included, fresh against the current
                # global bindings -- the same reason ppc_ensure_gramsuite_loaded!
                # itself uses `@eval import GRAMSuite` rather than a plain
                # `import`.
                load_gram && @eval ppc_ensure_gramsuite_loaded!()
                @eval ppc_startup_trace("loaded")
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
    # Resolve UNDER THE MODE'S ENVIRONMENT. OuterRouteTuning's defaults read
    # the profile switches from ENV (SPACEAGORA_PARALLEL_POLICY_V2 selects the
    # core-budget Monte Carlo default and the SPACEAGORA_PERF_PROCS worker cap),
    # and this used to run before ppc_run_sample_batch applied the mode env --
    # so every adaptive profile was routed with the shipped defaults whatever
    # its profile declared, and a paired R6-vs-R5 campaign probe measured two
    # copies of the same route. The shipped R4/R5 defaults read nothing from
    # ENV on the cold path, so their routing is unchanged by this.
    route = try withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=sample_count)...) do
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
        #
        # That freshness carries a second property, currently load-bearing and
        # easy to lose. An empty state means an empty snapshot, and
        # select_outer_route! guards its entire statistics section behind
        # `if !isempty(snapshot)`, so on a fresh state only the COLD path runs:
        # default_outer_route and the candidate filter. The
        # accumulated-statistics half -- the proven-default guard, the
        # confidence ranking, the split-width ladder, and the feedback recorder
        # this harness never calls at all -- is unreachable, and changes to it
        # cannot move a benchmark number. That is why work on those parts can
        # land mid-measurement-campaign without invalidating runs in flight.
        #
        # Be precise about which half, because the other half is NOT inert.
        # Changes to the cold path are fully visible to R4 and R5 and do move
        # published numbers: reversing the Monte Carlo default from process to
        # threads, and later sending it back to process at one thread, both
        # landed in default_outer_route and both moved every adaptive Monte
        # Carlo figure by tens of percent. "Routing changes cannot move a
        # benchmark number" is true of the statistics half and false of the
        # cold path, and a reader who takes the unqualified version will trust
        # a stale CSV.
        #
        # Persist route state across points to save time and that isolation is
        # gone: the harness would start measuring the guard's convergence rather
        # than the router's cold decision, points would stop being independent,
        # and their order would matter. If that is ever wanted it needs to be a
        # deliberate, separately reported configuration, not an optimisation.
        #
        # For a PINNED backend the guarantee is stronger still and does not rest
        # on the empty snapshot at all: the first line of this function is
        # `mode.backend == "auto" || return mode`, so a mode declaring
        # "threads", "process" or "serial" never reaches select_outer_route!.
        # The only other ParallelProfiles entry points the harness touches are
        # _machine_parallel_class, a hardware classifier holding no routing
        # state, and profile_env_pairs in modes.jl. So for R0 through R3 the
        # routing sources are not merely unreachable, they are never called --
        # which is what rules them out as the cause when a static route's timing
        # shifts between commits, and sends the search into the rest of src/
        # instead.
        PP.select_outer_route!(
            PP.OuterRouteState(),
            features;
            machine_class=PP._machine_parallel_class(),
            threads_available=Threads.nthreads() > 1,
            parallel_enabled=true,
        )
    end catch err
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
    # Adaptive profiles run their batch through the shipped campaign runner
    # (see ppc_run_adaptive_batch); the static modes below pin their dispatch.
    if sample_count > 1 && mode.backend == "auto" && mode.policy_adaptive && ppc_adaptive_via_runner()
        return ppc_run_adaptive_batch(case, cfg, mode, sample_indices, sample_seeds)
    end
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
    ppc_startup_trace("timed_start")
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
            # Literally what the campaign runner does after its own dispatch,
            # through the same function: collect on the idle workers between
            # batches so the next batch does not pay a mid-round stall. Calling
            # it rather than repeating it is what keeps the static process route
            # and the adaptive runner on the same footing, including under
            # SPACEAGORA_POOL_WORKER_GC.
            SpaceAGORA.SimulationCampaigns._collect_on_workers(worker_ids)
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
        # The mode as RESOLVED (backend auto -> threads/process), which is the
        # one the timed batch ran under. The caller records the row's env from
        # it; recording from the unresolved mode wrote backend=auto and hid
        # every backend-conditional pair (the inner budget among them).
        mode=mode,
        env_extra=Pair{String, String}[],
        policy=ppc_policy_columns(
            isempty(results) ? nothing :
                (results[end] isa NamedTuple && haskey(results[end], :policy) ? results[end].policy : nothing)
        )
    )
end

# SPACEAGORA_PPC_ADAPTIVE_VIA_RUNNER=0 restores the pinned dispatch for the
# adaptive modes (route resolved once, then this harness's own threads/process
# loop), for attribution runs against the runner.
@inline function ppc_adaptive_via_runner()::Bool
    return _ppc_bool(get(ENV, "SPACEAGORA_PPC_ADAPTIVE_VIA_RUNNER", "1"))
end

"""
    ppc_run_adaptive_batch(case, cfg, mode, sample_indices, sample_seeds)

Time an adaptive profile's sample batch through `run_monte_carlo(threads=:auto)`
-- the shipped campaign runner -- rather than this harness's own dispatch.

This harness used to resolve the adaptive route once (ppc_resolve_outer_backend)
and then run the batch itself: `Threads.@threads` for threads, `pmap` over its
pool for process. That measured the route DECISION but not the route's
EXECUTION, and the two diverged once V2 gave the process route mixed dispatch
(pool workers plus the coordinator's spare threads consuming one queue): a
harness `pmap` over W one-thread workers can only ever show W-way concurrency,
so B15's middle splits could not show the change at all. The runner is what a
user's campaign actually calls, so its dispatch is the thing to measure.

Preserved from the pinned path: the route decision is still cold (a fresh
OuterRouteState per batch), the same warm-ups run first, and the pool is
provisioned and warmed before the clock starts -- this harness's workers,
which carry the study files, are ADOPTED into the runner's pool rather than
letting it spawn a second set. The coordinator's mode env is applied with
outer_tasks=1, i.e. without declaring an outer split: the runner yields to an
enclosing split and would run serially, and it declares the split (and each
sample's inner budget) itself around the dispatch.

The sample function runs a coordinator-side sample directly and a worker-side
one through ppc_process_sample_task, whose own withenv is the worker's --
never a withenv inside a concurrent coordinator task, since ENV is
process-global.
"""
# One outer-route state per worker process, shared by a point's repeats.
#
# The adaptive policy is used as a sequence of campaigns of one shape on one
# machine, and that is how it is measured here: repeat 1 is its cold answer,
# the later repeats are what it learned from the earlier ones (a rounds tie
# spends repeat 2 on the other parallel arm -- select_outer_route!'s
# explore_tie -- and exploits from repeat 3). Points are separate processes
# and stay cold with respect to each other. The state path is pinned to a file
# no run writes (the runner only persists its own global state, never this
# one), so production history on the machine never leaks into a row.
const _PPC_ADAPTIVE_ROUTE_STATE = Ref{Any}(nothing)

function _ppc_adaptive_route_state()
    st = _PPC_ADAPTIVE_ROUTE_STATE[]
    if st === nothing
        st = SpaceAGORA.ParallelProfiles.OuterRouteState()
        _PPC_ADAPTIVE_ROUTE_STATE[] = st
    end
    return st
end

# The campaign closure, cached per point.
#
# `CachingPool` keys its worker-side cache on the IDENTITY of the function it
# is handed (see `_acquire_dispatch_runner` in
# src/simulation/campaigns/monte_carlo.jl), so a closure rebuilt for every
# repeat is a guaranteed cache miss on every worker, every time: the campaign's
# closure is serialized to each of them again and torn down again. Repeats of
# one point are the same campaign run again, and that is what they should cost.
#
# The key is what the sample path actually reads out of the config, not the
# config itself: `ppc_single_config` reads `profile`, `ppc_mode_env_pairs`
# reads `process_workers` and `solver_mode`, and `ppc_solve_once` reads
# nothing. The config object DOES differ between repeats of one point
# (`worker_seed`, `warmup`, `worker_repeat` -- see `ppc_run_worker_performance`),
# and the sample index and seed both arrive as the sample's own job rather than
# being derived from a captured `worker_seed`, so none of the three is read
# through the reused closure. The key is the statement that nothing the closure
# can act on differs; if a field the sample path reads is ever added to it, it
# belongs in the key.
#
# SPACEAGORA_PPC_SAMPLE_FN_REUSE=0 rebuilds the closure per batch, which is
# what this harness did before, for measuring the difference.
const _PPC_SAMPLE_FN_CACHE = Ref{Any}(nothing)

@inline function _ppc_sample_fn_reuse()::Bool
    return _ppc_bool(get(ENV, "SPACEAGORA_PPC_SAMPLE_FN_REUSE", "1"))
end

function _ppc_sample_fn(case_name::String, mode_name::String, cfg::PPCConfig)
    key = (case_name, mode_name, cfg.profile, cfg.process_workers, cfg.solver_mode)
    cached = _PPC_SAMPLE_FN_CACHE[]
    if _ppc_sample_fn_reuse() && cached !== nothing && cached.key == key
        return cached.fn
    end
    # One job is (index, seed): a coordinator-side sample and a worker-side one
    # need the index, and deriving it from a captured seed origin is what used
    # to tie this closure to one repeat.
    fn = job -> begin
        idx, seed = job
        Distributed.myid() == 1 ?
            ppc_run_sample_once(case_name, cfg, idx, seed) :
            ppc_process_sample_task(case_name, cfg, mode_name, idx, seed)
    end
    _PPC_SAMPLE_FN_CACHE[] = (key=key, fn=fn)
    return fn
end

function ppc_run_adaptive_batch(
    case::PPCCaseSpec, cfg::PPCConfig, mode::PPCModeSpec,
    sample_indices::Vector{Int}, sample_seeds::Vector{Int}
)
    sample_count = length(sample_seeds)
    SCamp = SpaceAGORA.SimulationCampaigns
    PP = SpaceAGORA.ParallelProfiles
    # Warm-ups as on the pinned path; the timed samples will see an active
    # outer split, and so should these.
    withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=sample_count)...) do
        for warmup_idx in 1:cfg.warmup
            warmup_args = ppc_single_config(case.name, cfg; seed=cfg.worker_seed - warmup_idx, mc_index=warmup_idx)
            try
                ppc_solve_once(warmup_args, cfg)
            catch
            end
        end
    end
    # Pool before the clock, adopted into the runner's pool. One worker can
    # never carry the process route (it withdraws below two), so none then.
    case_name = case.name
    worker_seed = cfg.worker_seed
    sample_jobs = [(sample_indices[i], sample_seeds[i]) for i in eachindex(sample_seeds)]
    sample_fn = _ppc_sample_fn(case_name, mode.name, cfg)
    if cfg.process_workers >= 2
        ids = ppc_ensure_process_workers!(cfg.process_workers)
        SpaceAGORA.adopt_process_workers!(SpaceAGORA.campaign_process_pool(), ids)
        # Warm each worker through the SAME path the runner will use -- the
        # sample closure wrapped by _run_monte_carlo_sample -- not just
        # ppc_process_sample_task: the closure's first call on a worker JITs
        # that wrapper, and measured at (2 workers, 6 threads) that first call
        # put 2.4 s on repeat 1 alone (5.45 s against 3.01 / 3.48 s).
        if cfg.warmup > 0 && !isempty(ids)
            @sync for (offset, w) in enumerate(ids)
                @async try
                    remotecall_wait(w, (offset, cfg.worker_seed - offset)) do warm_job
                        SpaceAGORA.SimulationCampaigns._run_monte_carlo_sample(sample_fn, 0, warm_job)
                    end
                catch
                end
            end
        end
    end
    # The probe config and route features are the route DECISION's inputs; the
    # pinned path resolves its route before the clock, so this does too.
    features = withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=1)...) do
        probe = ppc_single_config(case_name, cfg; seed=worker_seed, mc_index=1)
        SCamp.campaign_route_features(probe; samples=sample_count)
    end
    route_state = _ppc_adaptive_route_state()
    state_path = joinpath(tempdir(), "spaceagora_ppc_outer_route_state_$(getpid()).toml")
    # Warm the path that is about to be timed. The warm-ups above run single
    # solves on the coordinator and one sample per worker, which is the whole
    # of a pinned route's machinery but none of an adaptive one's: the campaign
    # planner, the dispatchers, the pool feeders and the result plumbing were
    # first compiled inside the first timed repeat, and the workers' heaps
    # first grew to a campaign's worth of garbage there too. Measured on TRX50
    # at P4's top budget, an adaptive route's first timed repeat read 5.3 s
    # against the pinned pool's 3.0 s after the same warm-ups, and its next
    # five carried the workers' catch-up collections. So each warm-up also
    # runs one untimed campaign through the runner, exactly as the timed
    # repeat will.
    for _ in 1:cfg.warmup
        withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=1)...,
                "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => state_path) do
            try
                SCamp.run_monte_carlo(sample_fn, sample_jobs; threads=:auto,
                                      route_features=features, route_state=route_state)
            catch
            end
        end
    end
    GC.gc()
    ppc_startup_trace("timed_start")
    batch_started = time()
    r = withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=1)...,
                "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => state_path) do
        SCamp.run_monte_carlo(sample_fn, sample_jobs; threads=:auto,
                              route_features=features, route_state=route_state)
    end
    batch_wall = Float64(time() - batch_started)
    results = Vector{Any}(undef, sample_count)
    for smp in r.samples
        # ppc_run_sample_once reports solver failures in its own row; a failure
        # surfacing HERE is infrastructure (a dead worker, a serialization
        # error) and the point should fail loudly, as the pinned path did.
        smp.success || throw(ErrorException(
            "adaptive batch sample $(smp.index) failed outside the solve: $(sprint(showerror, smp.error))"))
        results[smp.index] = smp.value
    end
    # What the runner declared to the coordinator-side samples, for the row.
    share = r.local_slots > 0 ? max(1, fld(Threads.nthreads(), r.local_slots)) :
            (r.route === :threads ? max(1, fld(Threads.nthreads(), max(1, r.threads))) : 0)
    env_extra = Pair{String, String}["SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"]
    share > 0 && push!(env_extra, "SPACEAGORA_INNER_THREAD_BUDGET" => string(share))
    scope = r.local_slots > 0 ?
        "adaptive_mixed_w$(r.threads - r.local_slots)_l$(r.local_slots)_sample_batch" :
        "adaptive_$(r.route)_sample_batch"
    return (
        results=results,
        batch_wall_time_s=batch_wall,
        actual_backend=string(r.route),
        execution_scope=scope,
        outer_tasks=r.threads,
        mode=mode,
        env_extra=env_extra,
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
        rhs_plan_scheduler="none",
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
        rhs_plan_scheduler=string(g(:rhs_plan_scheduler, "none")),
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
    env_string = ppc_effective_env_string(batch.mode, repeat_cfg; outer_tasks=batch.outer_tasks)
    if !isempty(batch.env_extra)
        env_string *= ";" * join(("$(k)=$(v)" for (k, v) in batch.env_extra), ";")
    end
    sample_results = batch.results
    total_success = all(r -> r.success, sample_results)
    sample_wall_sum = sum(r -> Float64(r.wall_time_s), sample_results)
    # GC time and allocation, summed over the batch's samples: the two numbers that
    # separate 'this route does more work' from 'this route allocates and the
    # collector stops every concurrent sample' -- same wall time, different cause.
    sample_gc_sum = sum(r -> Float64(get(r, :gc_time_s, 0.0)), sample_results; init=0.0)
    sample_alloc_mb = sum(r -> Float64(get(r, :alloc_bytes, 0)), sample_results; init=0.0) / 2^20
    final_result = sample_results[end]
    final_retcode = total_success ? string(final_result.retcode) : join(unique(string(r.retcode) for r in sample_results if !r.success), "|")
    # The per-sample rows carry the exception text; keep the distinct texts of the
    # failed samples on the aggregate row too, or a campaign that threw reports
    # only the exception TYPE and the message has to be re-run for. Bounded so a
    # long stack cannot blow up the CSV cell.
    final_error_message = total_success ? "" :
        first(join(unique(String(get(r, :error_message, "")) for r in sample_results
                          if !r.success && !isempty(get(r, :error_message, ""))), " | "), 2000)
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
        error_message=final_error_message,
        # The slowest sample and its share of the campaign. @timed's counters are
        # process-global, so without this the aggregate row cannot say whether
        # one sample hung or every sample ran slow -- the two explanations for
        # a 100x campaign that need opposite fixes.
        max_sample_wall_time_s=maximum(Float64(r.wall_time_s) for r in sample_results; init=0.0),
        max_sample_index=(isempty(sample_results) ? 0 :
            argmax([Float64(r.wall_time_s) for r in sample_results])),
        wall_time_s=batch.batch_wall_time_s,
        sample_wall_time_sum_s=sample_wall_sum,
        mean_sample_wall_time_s=sample_wall_sum / max(1, samples),
        sample_gc_time_sum_s=sample_gc_sum,
        sample_alloc_mb_sum=sample_alloc_mb,
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

# GC collector flags for the worker subprocess's own `julia` launch
# (`--gcthreads`, `--heap-size-hint`), read from env at worker-spawn time --
# NOT part of PPCConfig, so this stays a change to worker launch flags only,
# not to the config surface other parts of this study own.
#
# `--gcthreads` is opt-in and unset by default (`SPACEAGORA_PPC_WORKER_GCTHREADS`),
# reproducing the pre-existing argv exactly when not given. See
# docs/architecture/heap_contention.md ("Collector settings") for the grid
# this was measured against: on a 12-core/24-thread workstation, at 8
# threads, on the P3/P4 shapes' outer_threads route -- well below the
# 16+-thread regime where single-heap GC contention actually shows up; see
# the doc for what to run at wider thread counts.
#
# `--heap-size-hint` is ON by default (a derived value; see the block
# comment above `_default_pool_worker_heap_hint_bytes` in
# src/parallel/process/worker_pool.jl for the diagnostic and the formula).
# `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT=off` removes the flag entirely; any
# other value overrides the computed default and is passed through verbatim.
function _ppc_worker_gc_flags(cfg::PPCConfig)::Vector{String}
    flags = String[]
    gcthreads = get(ENV, "SPACEAGORA_PPC_WORKER_GCTHREADS", "")
    isempty(gcthreads) || push!(flags, "--gcthreads=$(gcthreads)")

    heap_hint = strip(get(ENV, "SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT", ""))
    if !isempty(heap_hint)
        lowercase(heap_hint) == "off" || push!(flags, "--heap-size-hint=$(heap_hint)")
    else
        # Default (env unset): the same derived hint the process pool uses
        # by default (see src/parallel/process/worker_pool.jl and
        # docs/architecture/heap_contention.md), reusing its exact formula
        # and constants rather than a second copy of them. The harness's own
        # per-point worker subprocess showed the identical unbounded-heap
        # pattern the pool diagnostic did: a 4096-sample serial point
        # (process_workers not in play) reached 46 GB by itself with nothing
        # bounding its collector. `cfg.process_workers` stands in for
        # "pool_size" here -- it is the concurrent-worker budget this run was
        # configured with, even for a point that itself runs serially.
        default_bytes = SpaceAGORA.ParallelProcess._default_pool_worker_heap_hint_bytes(max(cfg.process_workers, 1))
        push!(flags, "--heap-size-hint=$(SpaceAGORA.ParallelProcess._format_heap_size_hint(default_bytes))")
    end
    return flags
end

# ── Precompile workload ──────────────────────────────────────────────────────
#
# Every point is a fresh Julia process, and so is every process-pool worker, so
# each one pays package loading plus compilation of the solver specializations
# it is about to run before it can time anything. The paper harness's precompile
# workload (benchmarks/studies/paper_parallelization_benchmarks/workload, built
# with its build_workload.sh) holds those specializations in a package image;
# a worker that loads it skips most of that compilation. See the README there
# for what it covers and the measured effect.
#
# The workload's environment is STACKED behind the repository project rather
# than replacing it: a worker keeps --project=<repo>, so it loads the same
# SpaceAGORA image, resolves every other package the same way and has the same
# active project (which the runtime's process pool and state paths key on) as a
# worker launched without the workload. The environment only adds the workload
# package, found second on JULIA_LOAD_PATH.
#
#   SPACEAGORA_PPB_WORKLOAD      0/off (default): never use it; auto: use the
#                                image when it is built and current, otherwise
#                                warn and run without it; 1/on: fail if it is
#                                not built and current. Off by default because
#                                the image shifts timed repeats by a few percent
#                                (code placement; see the workload README).
#   SPACEAGORA_PPB_WORKLOAD_ENV  the workload environment (default
#                                output/paper_workload/env in this checkout).

ppc_workload_default_env() = joinpath(PPC_REPO_ROOT, "output", "paper_workload", "env")

# Cases a worker runs without the workload: the native-GRAM cases, which the
# workload leaves out (see ppb_workload_excluded in the workload's points.jl).
ppc_workload_skips_case(case::AbstractString) = any(c -> occursin(c, case), PPC_GRAM_LIVE_CASES)

# The worker's load path: the controller's, with the workload environment
# inserted right after the active project.
function ppc_workload_load_path(env::AbstractString)::String
    pathsep = Sys.iswindows() ? ";" : ":"
    raw = strip(get(ENV, "JULIA_LOAD_PATH", ""))
    entries = isempty(raw) ? ["@", "@v#.#", "@stdlib"] : String.(split(raw, pathsep))
    # A trailing empty entry in JULIA_LOAD_PATH means "append the defaults".
    if !isempty(entries) && isempty(last(entries))
        entries = vcat(entries[1:end-1], ["@", "@v#.#", "@stdlib"])
    end
    entries = unique(filter(!=(String(env)), entries))
    at = findfirst(==("@"), entries)
    at === nothing ? pushfirst!(entries, String(env)) : insert!(entries, at + 1, String(env))
    return join(entries, pathsep)
end

const _PPC_WORKLOAD_ENV_CHECKED = Dict{String, Bool}()

"""
    ppc_workload_env() -> Union{Nothing, String}

The workload environment the controller should launch workers with, or `nothing`
to launch them without it. Checked once per controller process and environment
path, in a subprocess started exactly the way a worker is (same project, same
stacked load path), by asking whether the workload package's image is current
for it -- which also fails when SpaceAGORA or any other dependency has changed
since the image was built.
"""
function ppc_workload_env()::Union{Nothing, String}
    setting = lowercase(strip(get(ENV, "SPACEAGORA_PPB_WORKLOAD", "0")))
    setting in ("0", "off", "false", "no") && return nothing
    required = setting in ("1", "on", "true", "yes")
    env_raw = strip(get(ENV, "SPACEAGORA_PPB_WORKLOAD_ENV", ""))
    env = abspath(isempty(env_raw) ? ppc_workload_default_env() : String(env_raw))
    fresh = get!(_PPC_WORKLOAD_ENV_CHECKED, env) do
        isfile(joinpath(env, "Project.toml")) || return false
        julia_bin = Base.julia_cmd().exec[1]
        probe = "exit(let id = Base.identify_package(\"$(PPC_WORKLOAD_PACKAGE)\"); " *
                "id !== nothing && Base.isprecompiled(id) end ? 0 : 3)"
        cmd = addenv(`$(julia_bin) --startup-file=no --project=$(PPC_REPO_ROOT) -e $(probe)`,
                     "JULIA_LOAD_PATH" => ppc_workload_load_path(env))
        return success(pipeline(cmd; stdout=devnull, stderr=devnull))
    end
    fresh && return env
    msg = "Precompile workload image not built or not current for $(env); " *
          "build it with benchmarks/studies/paper_parallelization_benchmarks/workload/build_workload.sh"
    required && throw(ErrorException(msg * " (SPACEAGORA_PPB_WORKLOAD=$(setting))."))
    setting == "auto" && isfile(joinpath(env, "Project.toml")) && @warn msg * "; running without it."
    return nothing
end

# Environment a worker for `case` is launched with when the workload is in use.
function ppc_workload_worker_env(workload_env::Union{Nothing, String}, case::AbstractString)
    (workload_env === nothing || ppc_workload_skips_case(case)) && return Pair{String, String}[]
    return [
        "JULIA_LOAD_PATH" => ppc_workload_load_path(workload_env),
        "SPACEAGORA_PPC_WORKLOAD" => "1",
        # The runtime's own process pool, should a campaign grow one here.
        "SPACEAGORA_PROCESS_WORKER_PRELOAD" => PPC_WORKLOAD_PACKAGE,
    ]
end

function ppc_worker_cmd(cfg::PPCConfig; case::String, mode::String, threads::Int, repeat::Int, seed::Int, mc_samples::Int, outfile::String, parity::Bool, repeats::Int=1,
                        workload_env::Union{Nothing, String}=nothing)
    julia_bin = Base.julia_cmd().exec[1]
    argv = String[
        julia_bin,
        "--threads=$(threads)",
        _ppc_worker_gc_flags(cfg)...,
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
    cmd = Cmd(_ppc_apply_cpu_pinning(argv, cfg.cpu_pinning, threads))
    extra = ppc_workload_worker_env(workload_env, case)
    return isempty(extra) ? cmd : addenv(cmd, extra...)
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
    # An empty parity list means "no parity", not "every case". ppc_resolve_cases
    # expands an empty request to the whole sorted catalog, which is what `cases`
    # wants -- ppc_parse_cli fills its profile default in first, so the controller
    # never sees an empty one -- but not what parity wants: a caller that builds a
    # PPCConfig directly has no other way to say "none". ppb_quick_phases does
    # exactly that, declaring parity_cases = String[] on all three of its phases,
    # and read literally a --quick run queued parity for every catalog case
    # against each non-serial mode. It stalled on atmo256_gram_live_10min, whose
    # parity point is two outer_tasks=1 solves of live native GRAM at N=256.
    parity_cases = isempty(cfg.parity_cases) ? String[] : ppc_resolve_cases(cfg.parity_cases)
    modes = ppc_mode_specs()
    unknown_modes = [m for m in cfg.modes if !haskey(modes, m)]
    isempty(unknown_modes) || throw(ArgumentError("Unknown mode(s): $(join(unknown_modes, ", "))"))

    # Refuse to time anything on a machine that is already loaded. See
    # ppc_assert_machine_quiet! -- this is the check whose absence let an entire
    # 115-point run be collected against another user's 26-hour job.
    ppc_assert_machine_quiet!()
    workload_env = ppc_workload_env()
    println("[parallelization-performance] profile=$(cfg.profile) outdir=$(outdir)")
    println("[parallelization-performance] load=$(round(ppc_load_average(); digits=2)) " *
            "headroom=$(round(100 * ppc_load_headroom(); digits=0))% " *
            "cores=$(_ppc_physical_core_count())")
    println("[parallelization-performance] cases=$(join(cases, ","))")
    println("[parallelization-performance] modes=$(join(cfg.modes, ","))")
    println("[parallelization-performance] precompile_workload=$(workload_env === nothing ? "off" : workload_env)")

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
                    parity=false,
                    workload_env=workload_env
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
                parity=true,
                workload_env=workload_env
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
