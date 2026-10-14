# Serial vs. 4-outer-process comparison for the mc_multisat Monte Carlo
# campaign (256 samples, 4 satellites/sample, 150km/600s LEO -- same
# per-sample shape as mc_multisat_thread_allocation.jl's known-good grid, just
# more samples), for both a gravity-only ("no_gram") and GRAM-surrogate
# ("surrogate") density model. N_SAMPLES is scaled up (not per-sample size)
# so the campaign is long enough (~seconds) that Distributed dispatch/IPC
# overhead doesn't swamp the real per-process parallel benefit -- an earlier
# attempt that instead scaled satellites/sample and mission length to do this
# hit a pathological non-convergent integration on some seed (17+ CPU-minutes
# on one sample) and was abandoned in favor of this safer axis.
#
# Follow-up to mc_multisat_thread_allocation.jl: that study found every
# outer_workers>=2 combo (Threads.@spawn-based) crashes for the GRAM
# surrogate/native models, root-caused to unlocked GRAM model *construction*
# racing across threads that share one process's address space (Finding 4,
# THREAD_ALLOCATION_AND_GRAM_CONCURRENCY_HANDOFF.md /
# project_gram_concurrent_construction_crash memory). Separate OS processes
# don't share that native global state, so this script routes outer
# parallelism through Distributed.addprocs instead of Threads.@spawn --
# the `outer_process` pattern already used in
# benchmarks/studies/parallelization_performance/execution.jl.
#
#   julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/mc_multisat_process_backend.jl
#
# Each Distributed worker is started with --threads=1 (no inner parallelism);
# this isolates the outer-process-vs-outer-thread comparison from inner
# thread-budget effects already covered by the earlier study.

using Distributed

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SAMPLE_LOGIC = joinpath(@__DIR__, "mc_multisat_process_sample.jl")
const MODES = ("no_gram", "surrogate")
const N_SATS_PER_SAMPLE = 4
const N_SAMPLES = 256
const OUTER_PROCESSES = 4
const OUT_CSV = joinpath(@__DIR__, "mc_multisat_process_backend_summary.csv")

include(SAMPLE_LOGIC)  # loads mc_run_sample/mc_build_sample_config on the driver too (serial route);
                        # also `using SpaceAGORA` transitively via examples/common.jl, so SpaceAGORA is
                        # a bare name here.

# Was a hand-rolled addprocs + serial-then-parallel bootstrap (worked around a
# precompile-cache pidfile race: N fresh processes hitting the
# SpaceAGORA+GRAMSuite+DataFrames+Plots extension-precompilation cascade at
# once reproducibly deadlocked one worker at 0% CPU on this machine/Julia
# 1.12.1). Now delegates process spawning + the SpaceAGORA/GRAMSuite/default
# SPICE bootstrap to SpaceAGORA.ParallelProcess's shared pool, which already
# bootstraps workers one at a time for the same reason -- see
# `parallel/process/worker_pool.jl`. Workers are shared with any other
# process-routed campaign in this session; only this study's own per-sample
# logic still needs including here.
const _MC_STUDY_BOOTSTRAPPED_WORKERS = Set{Int}()

function ensure_process_workers!(n::Int)
    pool = SpaceAGORA.campaign_process_pool()
    worker_ids = SpaceAGORA.ensure_process_workers!(pool, n)
    new_workers = [w for w in worker_ids if !(w in _MC_STUDY_BOOTSTRAPPED_WORKERS)]
    if !isempty(new_workers)
        for w in new_workers
            remotecall_wait(w, SAMPLE_LOGIC) do path
                Core.eval(Main, :(include($path)))
                nothing
            end
        end
        union!(_MC_STUDY_BOOTSTRAPPED_WORKERS, new_workers)
    end
    return worker_ids
end

function run_serial(mode::String)
    # `invokelatest` sidesteps a Julia 1.12 world-age error on first-ever GRAM
    # surrogate construction in this process: GRAMSuite lazily loads/defines
    # methods (e.g. `set_library!`) on first use, and this top-level driver
    # script's functions were already compiled against the pre-load world age.
    # Distributed workers don't need this -- each `pmap` dispatch re-resolves
    # against the latest world automatically.
    Base.invokelatest(mc_run_sample, mode, N_SATS_PER_SAMPLE, 0, 0)  # warmup (JIT), discarded
    GC.gc()
    start = time()
    samples = [Base.invokelatest(mc_run_sample, mode, N_SATS_PER_SAMPLE, i, i) for i in 1:N_SAMPLES]
    wall = time() - start
    return samples, wall
end

function run_outer_process(mode::String, n_workers::Int)
    worker_ids = ensure_process_workers!(n_workers)
    pool = CachingPool(worker_ids)
    pmap(pool, worker_ids) do _
        mc_run_sample(mode, N_SATS_PER_SAMPLE, 0, 0)  # per-worker warmup, discarded
    end
    start = time()
    samples = pmap(pool, 1:N_SAMPLES) do i
        mc_run_sample(mode, N_SATS_PER_SAMPLE, i, i)
    end
    wall = time() - start
    return samples, wall
end

function print_result(label::String, mode::String, samples, wall::Float64)
    n_ok = count(s -> s.success, samples)
    println("$(label) mode=$(mode): wall_time_s=$(round(wall; digits=4)) n_ok=$(n_ok)/$(length(samples))")
    for s in samples
        s.success || println("  sample index=$(s.index) seed=$(s.seed) FAILED: $(s.error)")
    end
end

function main()
    println("Serial vs $(OUTER_PROCESSES)-outer-process Monte Carlo comparison")
    println("N_SAMPLES=$(N_SAMPLES) sats/sample=$(N_SATS_PER_SAMPLE)")
    println()

    rows = String[]
    for mode in MODES
        println("=== mode=$(mode) : serial ===")
        s_samples, s_wall = run_serial(mode)
        print_result("serial", mode, s_samples, s_wall)
        push!(rows, "$(mode),serial,1,$(round(s_wall; digits=4))," *
            "$(count(s -> s.success, s_samples)),$(length(s_samples))")
        println()

        println("=== mode=$(mode) : outer_process (n=$(OUTER_PROCESSES)) ===")
        p_samples, p_wall = run_outer_process(mode, OUTER_PROCESSES)
        print_result("outer_process", mode, p_samples, p_wall)
        push!(rows, "$(mode),outer_process,$(OUTER_PROCESSES),$(round(p_wall; digits=4))," *
            "$(count(s -> s.success, p_samples)),$(length(p_samples))")
        println()

        speedup = s_wall / p_wall
        println("speedup (serial / outer_process): $(round(speedup; digits=3))x")
        println()
    end

    open(OUT_CSV, "w") do io
        println(io, "mode,backend,workers,wall_time_s,n_ok,n_samples")
        for r in rows
            println(io, r)
        end
    end
    println("Saved: $(OUT_CSV)")
end

main()
