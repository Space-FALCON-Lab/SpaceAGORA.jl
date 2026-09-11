# Monte Carlo thread-allocation study: 8 samples, each an independently jittered
# 4-satellite LEO constellation (150km, 600s mission), for a gravity-only
# ("no_gram"), GRAM-surrogate ("surrogate"), and native/real-GRAM ("standard")
# density model.
#
# Two outer-dispatch backends are compared:
#
#   threads -- run_monte_carlo(...; threads=N), Threads.@spawn-based outer dispatch
#              within one process. At each total OS thread count T in {1,2,4,8},
#              sweeps every (outer_workers, inner_threads) split with
#              outer_workers*inner_threads == T, to find which allocation of a fixed
#              thread budget between outer (across-sample) and inner (within-sample,
#              across-satellite/effector) parallelism performs best:
#                T=1: (1,1)
#                T=2: (2,1) outer-only, (1,2) inner-only
#                T=4: (4,1) outer-only, (2,2) balanced, (1,4) inner-only
#                T=8: (8,1) outer-only, (4,2), (2,4), (1,8) inner-only
#              "standard" (native GRAM) used to be excluded here entirely: it crashed
#              (SPICE state corruption) whenever 2+ samples ran concurrently on
#              different threads. Fixed this session (GRAM constructor lock +
#              GRAM/SPICE lock unification) -- included now to demonstrate that fix.
#              It should no longer crash, but native GRAM's own single process-wide
#              lock still caps it near ~1x speedup regardless of outer_workers.
#
#   process -- run_monte_carlo(...; threads=:auto, ...) routed to
#              SpaceAGORA.ParallelProcess's real OS-process backend (this session's
#              other new work). Sweeps outer_workers alone (PROCESS_WORKER_LADDER,
#              inner_threads fixed at 1 -- ParallelProcess workers are currently
#              --threads=1 with no per-worker inner-thread budget, and outer-only is
#              already known to be the best allocation for this workload shape, so
#              this isn't a real loss of coverage). Separate OS processes share no
#              native GRAM/CSPICE state, so "standard" should scale here even though
#              it's lock-limited under "threads".
#
#   julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/mc_multisat_thread_allocation.jl
#
# Saves full per-sample and per-combo summary results to CSV (see OUT_SUMMARY_CSV /
# OUT_SAMPLES_CSV below) so this sweep does not need to be rerun to analyze results.

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WORKER_SCRIPT = joinpath(@__DIR__, "mc_multisat_thread_allocation_worker.jl")
const MODES = ("no_gram", "surrogate", "standard")
const OUT_SUMMARY_CSV = joinpath(@__DIR__, "mc_multisat_thread_allocation_summary.csv")
const OUT_SAMPLES_CSV = joinpath(@__DIR__, "mc_multisat_thread_allocation_samples.csv")

function combos_for_budget(total_threads::Int)
    return [(outer, total_threads ÷ outer) for outer in 1:total_threads if total_threads % outer == 0]
end

const THREAD_BUDGETS = [1, 2, 4, 8]
const PROCESS_WORKER_LADDER = [1, 2, 4, 8]

struct SweepPoint
    mode::String
    outer_backend::String
    outer_workers::Int
    inner_threads::Int
    subprocess_threads::Int
end

function sweep_points()::Vector{SweepPoint}
    points = SweepPoint[]
    for mode in MODES
        for t in THREAD_BUDGETS, (outer, inner) in combos_for_budget(t)
            push!(points, SweepPoint(mode, "threads", outer, inner, t))
        end
        for outer in PROCESS_WORKER_LADDER
            # Coordinator only needs 1 thread: the process pool's own workers
            # (each --threads=1) provide the parallelism, not this process.
            push!(points, SweepPoint(mode, "process", outer, 1, 1))
        end
    end
    return points
end

struct SampleRow
    mode::String
    outer_backend::String
    outer_workers::Int
    inner_threads::Int
    sample_index::Int
    seed::Int
    success::Bool
    elapsed_s::Float64
end

struct SummaryRow
    mode::String
    outer_backend::String
    outer_workers::Int
    inner_threads::Int
    wall_time_s::Union{Nothing, Float64}
    n_ok::Union{Nothing, Int}
    n_samples::Union{Nothing, Int}
    output::String
end

# Some (outer_workers>1, inner_threads>1) combos hit a known, pre-existing nested
# outer+inner thread-contention pathology (see this folder's handoff doc) that isn't
# a bounded slowdown -- e.g. mode=no_gram outer=4 inner=2 was observed to still be
# running after 16+ minutes of wall time for what takes 0.036s at outer=4 inner=1,
# far past even that doc's already-extreme 155x-1400x characterization. Without a
# per-point ceiling, one such combo can block the entire sweep indefinitely. A killed
# point is recorded as a timeout (wall_time_s=nothing, same as any other
# summary-parse failure) rather than left to hang.
const POINT_TIMEOUT_S = parse(Float64, get(ENV, "SPACEAGORA_MCTA_POINT_TIMEOUT_S", "300"))

function run_worker(p::SweepPoint)::Tuple{SummaryRow, Vector{SampleRow}}
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([
        julia_bin, "--project=$(REPO_ROOT)", "--threads=$(p.subprocess_threads)",
        WORKER_SCRIPT, p.mode, string(p.outer_workers), string(p.inner_threads), p.outer_backend,
    ])
    io = IOBuffer()
    proc = run(pipeline(cmd; stdout=io, stderr=io); wait=false)
    # Poll process_running + a wall-clock deadline instead of `wait(proc)`: observed
    # in practice that `wait(proc)` can stay blocked well past a killed subprocess's
    # own exit if it left orphaned children (e.g. SpaceAGORA.ParallelProcess pool
    # workers from an outer_backend="process" point) holding the write end of this
    # `pipeline(...)`'s stdout/stderr pipe open -- the OS-level pipe has no EOF to
    # give until every holder of the write end closes it, and `wait`/the pipe-drain
    # machinery underneath `run(pipeline(...))` doesn't resolve until then. An
    # explicit poll loop only depends on process_running, not pipe EOF.
    deadline = time() + POINT_TIMEOUT_S
    timed_out = false
    while process_running(proc)
        if time() >= deadline
            timed_out = true
            # Kill this point's own pool workers (direct children of proc) BEFORE
            # killing proc itself -- if proc dies first, its children are the ones
            # left holding the pipe open with no owner left to reap them via the
            # worker script's normal shutdown_process_pool! exit path (which a kill
            # bypasses entirely).
            try
                run(`pkill -9 -P $(proc.pid)`)
            catch
            end
            kill(proc, Base.SIGKILL)
            break
        end
        sleep(0.5)
    end
    # `wait(proc)` used to guarantee the pipeline's async stdout/stderr-draining task
    # had fully caught up before returning; polling process_running doesn't. A short
    # grace period avoids losing the last buffered chunk (e.g. the SUMMARY line) to
    # that race on the normal, non-timeout exit path.
    sleep(0.3)
    output = String(take!(io))
    if timed_out
        output *= "\nTIMEOUT: point exceeded $(POINT_TIMEOUT_S)s wall-clock budget and was killed " *
                   "(likely the nested outer+inner thread-contention pathology -- see README)."
    end

    samples = SampleRow[]
    for m in eachmatch(r"sample index=(\d+) seed=(\d+) success=(true|false) elapsed_s=([\d.]+)", output)
        push!(samples, SampleRow(
            p.mode, p.outer_backend, p.outer_workers, p.inner_threads,
            parse(Int, m.captures[1]), parse(Int, m.captures[2]),
            m.captures[3] == "true", parse(Float64, m.captures[4]),
        ))
    end

    sm = match(r"SUMMARY .*?wall_time_s=([\d.]+) n_ok=(\d+)/(\d+)", output)
    summary = if sm === nothing
        SummaryRow(p.mode, p.outer_backend, p.outer_workers, p.inner_threads, nothing, nothing, nothing, output)
    else
        SummaryRow(
            p.mode, p.outer_backend, p.outer_workers, p.inner_threads,
            parse(Float64, sm.captures[1]), parse(Int, sm.captures[2]), parse(Int, sm.captures[3]), output,
        )
    end
    return summary, samples
end

function run_sweep()
    summaries = SummaryRow[]
    all_samples = SampleRow[]
    points = sweep_points()
    total = length(points)
    done = 0
    for p in points
        done += 1
        print("[$(done)/$(total)] mode=$(p.mode) outer_backend=$(p.outer_backend) outer_workers=$(p.outer_workers) inner_threads=$(p.inner_threads) ... ")
        flush(stdout)
        summary, samples = run_worker(p)
        push!(summaries, summary)
        append!(all_samples, samples)
        if summary.wall_time_s === nothing
            println("FAILED to parse summary; raw output follows:")
            println(summary.output)
        else
            println("$(summary.wall_time_s) s  ($(summary.n_ok)/$(summary.n_samples) ok)")
        end
    end
    return summaries, all_samples
end

function save_csv(path::String, header::String, rows::Vector{String})
    open(path, "w") do io
        println(io, header)
        for r in rows
            println(io, r)
        end
    end
    println("Saved: $(path)")
end

function print_summary_table(summaries::Vector{SummaryRow})
    println()
    for mode in MODES
        println("mode=$(mode)")
        println(
            rpad("outer_backend", 15), rpad("outer", 7), rpad("inner", 7),
            rpad("wall_time_s", 13), rpad("speedup", 10), "n_ok"
        )
        # Shared baseline per mode: the "threads" backend's single-worker point
        # (a plain in-process serial run). Reused for both backends' speedup
        # columns so "process" isn't graded against its own pool-spinup overhead.
        baseline_row = only(filter(
            s -> s.mode == mode && s.outer_backend == "threads" && s.outer_workers == 1 && s.inner_threads == 1,
            summaries
        ))
        baseline = baseline_row.wall_time_s
        for backend in ("threads", "process")
            for s in filter(s -> s.mode == mode && s.outer_backend == backend, summaries)
                speedup = (baseline === nothing || s.wall_time_s === nothing) ? NaN : baseline / s.wall_time_s
                println(
                    rpad(s.outer_backend, 15), rpad(string(s.outer_workers), 7), rpad(string(s.inner_threads), 7),
                    rpad(s.wall_time_s === nothing ? "FAILED" : string(round(s.wall_time_s; digits=4)), 13),
                    rpad(isnan(speedup) ? "-" : string(round(speedup; digits=3)), 10),
                    s.n_ok === nothing ? "-" : "$(s.n_ok)/$(s.n_samples)",
                )
            end
        end
        println()
    end
end

function main()
    points = sweep_points()
    println("Monte Carlo thread-allocation sweep: modes=$(MODES)")
    println("threads backend: thread budgets=$(THREAD_BUDGETS); process backend: outer-worker ladder=$(PROCESS_WORKER_LADDER)")
    println("Total sweep points: $(length(points)) subprocess launches (process-backend points additionally pay one-time worker-process precompilation on their first run).")
    println()

    summaries, samples = run_sweep()
    print_summary_table(summaries)

    save_csv(
        OUT_SUMMARY_CSV,
        "mode,outer_backend,outer_workers,inner_threads,wall_time_s,n_ok,n_samples",
        [
            "$(s.mode),$(s.outer_backend),$(s.outer_workers),$(s.inner_threads)," *
            "$(s.wall_time_s === nothing ? "" : s.wall_time_s),$(s.n_ok === nothing ? "" : s.n_ok)," *
            "$(s.n_samples === nothing ? "" : s.n_samples)"
            for s in summaries
        ],
    )
    save_csv(
        OUT_SAMPLES_CSV,
        "mode,outer_backend,outer_workers,inner_threads,sample_index,seed,success,elapsed_s",
        [
            "$(r.mode),$(r.outer_backend),$(r.outer_workers),$(r.inner_threads)," *
            "$(r.sample_index),$(r.seed),$(r.success),$(r.elapsed_s)"
            for r in samples
        ],
    )
end

main()
