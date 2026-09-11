# Entry point for mc_multisat_thread_allocation.jl -- runs ONE (mode, outer_backend,
# outer_workers, inner_threads) configuration: an 8-sample Monte Carlo campaign where
# each sample is an independently jittered 4-satellite LEO constellation (150km, 600s
# mission), propagated as one coupled run_simulation call via SpaceAGORA.run_monte_carlo.
# Sample-building logic (build_sample_config/env_pairs) lives in the companion
# mc_multisat_thread_allocation_sample.jl, which this file includes and which is also
# `include`d directly on SpaceAGORA.ParallelProcess pool workers for outer_backend="process"
# (see ensure_process_workers_with_sample_logic! below) -- kept ARGS-independent for
# exactly that reason (a pool worker's own ARGS is empty).
#
#   julia --project=. --threads=<T> mc_multisat_thread_allocation_worker.jl \
#       <mode> <outer_workers> <inner_threads> <outer_backend>
#
# mode is one of:
#   no_gram   -- NoAtmosphereModel, gravity (L20 harmonics) only, no shared native
#                library/lock -- expected to scale well.
#   surrogate -- GRAM offline surrogate (precomputed interpolation, no native calls,
#                no lock) -- the "GRAM case" for this study.
#   standard  -- native/real GRAM point density. Previously excluded from this study
#                entirely: it crashed (SPICE state corruption) whenever 2+ independent
#                run_simulation calls executed concurrently on different threads,
#                exactly what outer_backend="threads" with outer_workers>1 does. Fixed
#                this session (GRAM constructor lock + GRAM/SPICE lock unification) --
#                included here specifically to demonstrate that fix, and to show
#                outer_backend="process" (a separate OS process per sample -- no
#                shared native GRAM/CSPICE state to race on in the first place) scaling
#                past native GRAM's own single process-wide lock, which still caps
#                outer_backend="threads" near ~1x regardless of outer_workers.
#
# outer_backend is one of:
#   threads -- run_monte_carlo(...; threads=outer_workers), the original fixed-count
#              API. T (--threads=<T>) must be >= outer_workers*inner_threads.
#   process -- run_monte_carlo(...; threads=:auto, ...) routed to
#              SpaceAGORA.ParallelProcess's real OS-process backend via a route_tuning
#              that pins process_max_workers=outer_workers and forces the :process
#              route, with a fresh OuterRouteState() so no bandit exploration history
#              carries over between sweep points. inner_threads is fixed at 1 by the
#              controller for this backend (ParallelProcess workers are currently
#              --threads=1 with no per-worker inner-thread budget), and T=1 is
#              sufficient -- the process pool's own workers provide the parallelism,
#              not this coordinator process.
#
# outer_workers is the outer worker-task/process count (how many of the 8 samples run
# concurrently). inner_threads is the thread budget available to each sample's OWN
# within-sample parallelism (across the 4 satellites/effectors of that one sample),
# via SPACEAGORA_INNER_THREAD_BUDGET -- meaningful for outer_backend="threads" only.

using Distributed

include(joinpath(@__DIR__, "mc_multisat_thread_allocation_sample.jl"))

mode = length(ARGS) >= 1 ? ARGS[1] : ""
mode in ("no_gram", "surrogate", "standard") || error(
    "Usage: julia --project=. --threads=<N> $(basename(@__FILE__)) <no_gram|surrogate|standard> <outer_workers> <inner_threads> <threads|process>"
)
const OUTER_WORKERS = parse(Int, ARGS[2])
const INNER_THREADS = parse(Int, ARGS[3])
const OUTER_BACKEND = length(ARGS) >= 4 ? ARGS[4] : "threads"
OUTER_BACKEND in ("threads", "process") || error("outer_backend must be \"threads\" or \"process\", got \"$(OUTER_BACKEND)\"")

# Bootstraps SpaceAGORA.ParallelProcess pool workers with THIS study's sample-building
# logic (build_sample_config/env_pairs) -- the pool's generic bootstrap (using
# SpaceAGORA + GRAMSuite + default SPICE kernels) has no way to know about it. Mirrors
# mc_multisat_process_backend.jl's ensure_process_workers!. Pre-warm before dispatch so
# run_monte_carlo's own internal ensure_process_workers! call (inside the adaptive
# :auto path) just finds already-bootstrapped workers.
function ensure_process_workers_with_sample_logic!(n::Int)::Vector{Int}
    pool = SpaceAGORA.campaign_process_pool()
    worker_ids = SpaceAGORA.ensure_process_workers!(pool, n)
    sample_file = joinpath(@__DIR__, "mc_multisat_thread_allocation_sample.jl")
    for w in worker_ids
        remotecall_wait(w, sample_file) do path
            Core.eval(Main, :(include($path)))
            nothing
        end
    end
    return worker_ids
end

function run_campaign()
    pairs = env_pairs(mode, OUTER_BACKEND, OUTER_WORKERS, INNER_THREADS)
    return withenv(pairs...) do
        if OUTER_BACKEND == "threads"
            SpaceAGORA.run_monte_carlo(1:N_SAMPLES; threads=OUTER_WORKERS) do seed
                args = build_sample_config(mode, seed)
                result = SpaceAGORA.run_simulation(
                    args; isolate_state=false, return_solution=true, return_solver_metadata=true
                )
                result.solution
            end
        else
            ensure_process_workers_with_sample_logic!(OUTER_WORKERS)
            tuning = SpaceAGORA.OuterRouteTuning(
                mc_process_min_samples=1,
                mc_process_min_mission_s=0.0,
                process_max_workers=OUTER_WORKERS,
            )
            SpaceAGORA.run_monte_carlo(
                1:N_SAMPLES;
                threads=:auto,
                route_features=SpaceAGORA.campaign_route_features(
                    samples=N_SAMPLES, n_sats=N_SATS_PER_SAMPLE, density_family=mode, mission_time_s=MISSION_TIME_S
                ),
                route_state=SpaceAGORA.OuterRouteState(),
                route_tuning=tuning
            ) do seed
                args = build_sample_config(mode, seed)
                result = SpaceAGORA.run_simulation(
                    args; isolate_state=false, return_solution=true, return_solver_metadata=true
                )
                result.solution
            end
        end
    end
end

function main()
    # Warmup (JIT) -- discarded, not timed.
    warmup = run_campaign()
    n_ok_warmup = length(warmup.successful)
    println("warmup: $(n_ok_warmup)/$(N_SAMPLES) samples succeeded")
    flush(stdout)
    n_ok_warmup == N_SAMPLES || @warn "warmup: $(length(warmup.failed)) sample(s) failed -- see raw output above."

    GC.gc()
    result = run_campaign()
    n_ok = length(result.successful)
    println("campaign: $(n_ok)/$(N_SAMPLES) samples succeeded, wall_time_s=$(round(result.elapsed_s; digits=4))")
    for s in result.samples
        println("sample index=$(s.index) seed=$(s.seed) success=$(s.success) elapsed_s=$(round(s.elapsed_s; digits=4))")
    end
    println(
        "SUMMARY mode=$(mode) outer_backend=$(OUTER_BACKEND) outer_workers=$(OUTER_WORKERS) inner_threads=$(INNER_THREADS) " *
        "total_threads=$(Threads.nthreads()) wall_time_s=$(round(result.elapsed_s; digits=4)) n_ok=$(n_ok)/$(N_SAMPLES)"
    )
    # This process is a one-shot subprocess launched per sweep point (see the
    # controller's run_worker) -- its own exit doesn't automatically clean up any
    # SpaceAGORA.ParallelProcess pool workers it spawned via ensure_process_workers!
    # (they're plain addprocs'd OS processes, not children tied to this process's own
    # lifetime in a way Julia/the OS reaps automatically). Left alone, each
    # outer_backend="process" point leaks its pool workers, and -- observed in
    # practice -- an orphaned pool worker can keep this process's own inherited
    # stdout/stderr pipe open after this process exits, hanging the controller's
    # `wait` on it indefinitely. Shut the pool down explicitly before exiting.
    if OUTER_BACKEND == "process"
        SpaceAGORA.shutdown_process_pool!(SpaceAGORA.campaign_process_pool())
    end
    return nothing
end

main()
