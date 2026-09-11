# Outer-backend (threads vs. process) comparison for run_constellation_ensemble at a
# FIXED LEO constellation size, native GRAM only -- the workload class fixed by this
# session's GRAM-constructor-lock and GRAM/SPICE-lock-unification fixes. Previously:
# crashed outright at ensemble outer_workers>=2 for native GRAM. Even once that's
# fixed, native GRAM's own single process-wide lock caps thread-based scaling near
# ~1x regardless of worker count (see leo_thread_scaling_by_mode.jl's own "standard"
# points). SpaceAGORA.ParallelProcess's real OS-process backend (this session's other
# new work) sidesteps that lock entirely: separate OS processes don't share it.
#
#   julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_ensemble_outer_backend_scaling.jl
#
# Override constellation size / outer-worker ladder:
#   SPACEAGORA_EOB_N_SATS=32 SPACEAGORA_EOB_WORKERS=1,2,4,8 julia --project=. \
#       benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_ensemble_outer_backend_scaling.jl
#
# Each (backend, outer_workers) point is its own subprocess, matching
# leo_thread_scaling_by_mode.jl's pattern (Julia's thread count is fixed at process
# startup, so the "threads" backend needs a fresh process per point regardless).
# "process" points additionally spin up their own nested worker processes; those
# workers' one-time SpaceAGORA/GRAMSuite precompilation is cached machine-wide
# (~/.julia/compiled), so only the first "process" point in a cold run pays the full
# cost -- later points and later runs reuse the compiled cache.

using Plots

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WORKER_SCRIPT = joinpath(@__DIR__, "leo_ensemble_outer_backend_scaling_worker.jl")
const N_SATS = parse(Int, get(ENV, "SPACEAGORA_EOB_N_SATS", "64"))
const WORKER_LADDER = [parse(Int, strip(tok)) for tok in split(get(ENV, "SPACEAGORA_EOB_WORKERS", "1,2,4,8"), ",")]
const BACKENDS = ("threads", "process")
const BACKEND_LABELS = Dict(
    "threads" => "Threads.@spawn (in-process, lock-limited)",
    "process" => "OS processes (SpaceAGORA.ParallelProcess)",
)

struct BackendResult
    backend::String
    outer_workers::Int
    wall_s::Union{Nothing, Float64}
    n_ok::Union{Nothing, Int}
    output::String
end

function run_worker(backend::String, outer_workers::Int)::BackendResult
    julia_bin = Base.julia_cmd().exec[1]
    # The "threads" backend needs a matching Julia thread pool at process startup.
    # The "process" backend's real parallelism lives in its own worker processes
    # (each --threads=1, spun up by SpaceAGORA.ParallelProcess), so this coordinator
    # process only ever needs 1.
    coordinator_threads = backend == "threads" ? outer_workers : 1
    cmd = Cmd([
        julia_bin, "--project=$(REPO_ROOT)", "--threads=$(coordinator_threads)",
        WORKER_SCRIPT, string(N_SATS), string(outer_workers), backend,
    ])
    io = IOBuffer()
    try
        run(pipeline(cmd; stdout=io, stderr=io); wait=true)
    catch e
        e isa ProcessFailedException || rethrow()
    end
    output = String(take!(io))
    m = match(r"median wall time \(.*?\):\s*([\d.]+)\s*s", output)
    wall_s = m === nothing ? nothing : parse(Float64, m.captures[1])
    ok_m = match(r"repeat 1 \([\d.]+ s\): (\d+)/(\d+) members succeeded", output)
    n_ok = ok_m === nothing ? nothing : parse(Int, ok_m.captures[1])
    return BackendResult(backend, outer_workers, wall_s, n_ok, output)
end

function run_sweep()::Vector{BackendResult}
    results = BackendResult[]
    total = length(BACKENDS) * length(WORKER_LADDER)
    done = 0
    for backend in BACKENDS, outer_workers in WORKER_LADDER
        done += 1
        print("[$(done)/$(total)] backend=$(backend) outer_workers=$(outer_workers) ... ")
        flush(stdout)
        r = run_worker(backend, outer_workers)
        push!(results, r)
        if r.wall_s === nothing
            println("FAILED to parse wall time; raw output follows:")
            println(r.output)
        else
            println("$(round(r.wall_s; digits=4)) s (n_ok=$(r.n_ok)/$(N_SATS))")
        end
    end
    return results
end

function series_for(results::Vector{BackendResult}, backend::String)
    xs = Int[]
    ys = Float64[]
    for outer_workers in WORKER_LADDER
        r = only(filter(r -> r.backend == backend && r.outer_workers == outer_workers, results))
        r.wall_s === nothing && continue
        push!(xs, outer_workers)
        push!(ys, r.wall_s)
    end
    return xs, ys
end

function print_summary(results::Vector{BackendResult})
    println()
    println("Summary (N_SATS=$(N_SATS), native GRAM, ensemble route, outer_workers=$(WORKER_LADDER)):")
    for backend in BACKENDS
        xs, ys = series_for(results, backend)
        isempty(ys) && continue
        baseline = ys[1]
        println()
        println("backend=$(backend) ($(BACKEND_LABELS[backend]))")
        println(rpad("outer_workers", 15), rpad("wall_s", 12), rpad("speedup", 10), "efficiency")
        for (w, y) in zip(xs, ys)
            speedup = baseline / y
            efficiency = speedup / w
            println(
                rpad(string(w), 15),
                rpad(string(round(y; digits=4)), 12),
                rpad(string(round(speedup; digits=3)), 10),
                string(round(efficiency; digits=3)),
            )
        end
    end
end

function save_csv(results::Vector{BackendResult})
    path = joinpath(@__DIR__, "leo_ensemble_outer_backend_scaling_summary.csv")
    open(path, "w") do io
        println(io, "backend,outer_workers,wall_s,n_ok,n_sats")
        for r in results
            println(io, "$(r.backend),$(r.outer_workers),$(r.wall_s === nothing ? "" : r.wall_s)," *
                "$(r.n_ok === nothing ? "" : r.n_ok),$(N_SATS)")
        end
    end
    println("Saved: $(path)")
end

function make_plot(results::Vector{BackendResult})
    plt = plot(
        xlabel="Outer workers",
        ylabel="Wall time (s)",
        title="Ensemble outer-backend scaling: threads vs. process\n(LEO, $(N_SATS) sats, native GRAM, 600s mission)",
        xscale=:log2,
        yscale=:log10,
        legend=:topright,
        size=(750, 550),
        margin=5Plots.mm,
    )
    markers = Dict("threads" => :diamond, "process" => :circle)
    for backend in BACKENDS
        xs, ys = series_for(results, backend)
        isempty(ys) && continue
        plot!(plt, xs, ys; label=BACKEND_LABELS[backend], marker=markers[backend])
    end
    out_path = joinpath(@__DIR__, "leo_ensemble_outer_backend_scaling.png")
    savefig(plt, out_path)
    println()
    println("Saved: $(out_path)")
    return plt
end

function main()
    println("Ensemble outer-backend scaling sweep: N_SATS=$(N_SATS), backends=$(BACKENDS), outer_worker ladder=$(WORKER_LADDER)")
    println("Estimated wall time: several minutes (process-backend points additionally pay one-time worker-process precompilation on their first run).")
    println()

    results = run_sweep()
    print_summary(results)
    save_csv(results)
    make_plot(results)
end

main()
