# Thread-count scaling sweep for a fixed-size LEO constellation, run separately
# for each of the three density-model modes supported by
# leo_constellation_size_scaling_worker.jl:
#
#   no_gram   -- NoAtmosphereModel, gravity-only (no shared native library, no lock)
#   surrogate -- GRAM offline surrogate (precomputed interpolation, no native calls, no lock)
#   standard  -- real/native GRAM, no lookahead cache (single global lock inside the
#                vendored native library -- see LEO_GRAM_CONSTELLATION_HANDOFF.md)
#
# Purpose: isolate how much of the "threads don't help much" finding from the
# 2048-sat/1hr native-GRAM benchmark (LEO_GRAM_CONSTELLATION_HANDOFF.md) is
# specific to native GRAM's global lock vs. a general property of constellation
# parallelism in this codebase. Julia's thread count is fixed at process
# startup, so each (mode, threads) point is its own subprocess via the existing
# leo_constellation_size_scaling_worker.jl.
#
#   julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_thread_scaling_by_mode.jl
#
# Override satellite count / thread ladder:
#   SPACEAGORA_TS_N_SATS=512 SPACEAGORA_TS_THREADS=1,2,4,8,16 julia --project=. \
#       benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_thread_scaling_by_mode.jl

using Plots

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WORKER_SCRIPT = joinpath(@__DIR__, "leo_constellation_size_scaling_worker.jl")
const N_SATS = parse(Int, get(ENV, "SPACEAGORA_TS_N_SATS", "512"))
const THREAD_LADDER = [parse(Int, strip(tok)) for tok in split(get(ENV, "SPACEAGORA_TS_THREADS", "1,2,4,8"), ",")]
# ROUTE selects the execution route tested at each thread count:
#   monolithic -- one coupled state vector for all N_SATS satellites (flat RHS queue).
#   ensemble   -- run_constellation_ensemble: each satellite is an independent
#                 single-satellite run dispatched to an outer worker task.
const ROUTE = get(ENV, "SPACEAGORA_TS_ROUTE", "monolithic")
ROUTE in ("monolithic", "ensemble") || error("SPACEAGORA_TS_ROUTE must be \"monolithic\" or \"ensemble\", got \"$(ROUTE)\"")
# "surrogate_on" reuses mode="surrogate" but forces the density-callback thread gate
# open (see leo_constellation_size_scaling_worker.jl) instead of leaving it on "auto",
# to test whether the lock-free surrogate scales better once that gate isn't in the way.
# Only meaningful under the monolithic route -- the ensemble route always runs each
# (single-satellite) member with density-callback threading off, since there is no
# per-member benefit to gate in the first place.
const MODES = ROUTE == "monolithic" ?
    ("no_gram", "surrogate", "surrogate_on", "standard") : ("no_gram", "surrogate", "standard")
const MODE_LABELS = Dict(
    "no_gram" => "gravity-only (no atmosphere)",
    "surrogate" => "GRAM surrogate (density callback: auto)",
    "surrogate_on" => "GRAM surrogate (density callback: forced on)",
    "standard" => "native GRAM (no lookahead)",
)
const WORKER_MODE = Dict(
    "no_gram" => "no_gram", "surrogate" => "surrogate", "surrogate_on" => "surrogate", "standard" => "standard",
)

struct ThreadResult
    mode::String
    threads::Int
    median_s::Union{Nothing, Float64}
    output::String
end

function run_worker(mode::String, threads::Int)::ThreadResult
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([
        julia_bin, "--project=$(REPO_ROOT)", "--threads=$(threads)",
        WORKER_SCRIPT, string(N_SATS), WORKER_MODE[mode], ROUTE,
    ])
    env_overrides = mode == "surrogate_on" ?
        Dict("SPACEAGORA_TS_SURROGATE_DENSITY_CALLBACK_OVERRIDE" => "on") : Dict{String, String}()
    io = IOBuffer()
    # The ensemble route can crash a worker outright (e.g. native GRAM under 2+
    # concurrent ensemble members -- a shared-native-state hazard distinct from the
    # monolithic route's density-callback gate, and not covered by this worker's
    # existing FREEZE_PER_STEP/VACUUM_GRAM_CACHE=0 safeguards, which only protect
    # concurrency *within* one coupled solve). Catch that instead of aborting the
    # whole sweep so the rest of the ladder/modes still complete.
    try
        withenv(env_overrides...) do
            run(pipeline(cmd; stdout=io, stderr=io); wait=true)
        end
    catch e
        e isa ProcessFailedException || rethrow()
    end
    output = String(take!(io))
    m = match(r"median wall time \(.*?\):\s*([\d.]+)\s*s", output)
    median_s = m === nothing ? nothing : parse(Float64, m.captures[1])
    return ThreadResult(mode, threads, median_s, output)
end

function run_sweep()::Vector{ThreadResult}
    results = ThreadResult[]
    total = length(MODES) * length(THREAD_LADDER)
    done = 0
    for mode in MODES, threads in THREAD_LADDER
        done += 1
        print("[$(done)/$(total)] mode=$(mode) threads=$(threads) ... ")
        flush(stdout)
        r = run_worker(mode, threads)
        push!(results, r)
        if r.median_s === nothing
            println("FAILED to parse median wall time; raw output follows:")
            println(r.output)
        else
            println("$(round(r.median_s; digits=4)) s")
        end
    end
    return results
end

function series_for(results::Vector{ThreadResult}, mode::String)
    xs = Int[]
    ys = Float64[]
    for threads in THREAD_LADDER
        r = only(filter(r -> r.mode == mode && r.threads == threads, results))
        r.median_s === nothing && continue
        push!(xs, threads)
        push!(ys, r.median_s)
    end
    return xs, ys
end

function print_summary(results::Vector{ThreadResult})
    println()
    println("Summary (route=$(ROUTE), N_SATS=$(N_SATS), threads=$(THREAD_LADDER)):")
    for mode in MODES
        xs, ys = series_for(results, mode)
        isempty(ys) && continue
        baseline = ys[1]
        println()
        println("mode=$(mode) ($(MODE_LABELS[mode]))")
        println(rpad("threads", 10), rpad("median_s", 12), rpad("speedup", 10), "efficiency")
        for (t, y) in zip(xs, ys)
            speedup = baseline / y
            efficiency = speedup / t
            println(
                rpad(string(t), 10),
                rpad(string(round(y; digits=4)), 12),
                rpad(string(round(speedup; digits=3)), 10),
                string(round(efficiency; digits=3)),
            )
        end
    end
end

function make_plot(results::Vector{ThreadResult})
    plt = plot(
        xlabel="Threads",
        ylabel="Median wall time (s)",
        title="Thread-count scaling by density-model mode ($(ROUTE))\n(LEO, $(N_SATS) sats, 600s mission)",
        xscale=:log2,
        yscale=:log10,
        # "monolithic" has all four series clustered high (auto/on-surrogate and
        # native GRAM sit at 5-15s across the whole thread range) with only
        # gravity-only running low, leaving a trace-free band in the middle of the
        # plot -- :right lands the legend there. "ensemble" has native GRAM as a
        # single high-left point with everything else low, so :topright is empty.
        legend=(ROUTE == "monolithic" ? :right : :topright),
        size=(750, 550),
        margin=5Plots.mm,
    )
    markers = Dict("no_gram" => :circle, "surrogate" => :square, "surrogate_on" => :utriangle, "standard" => :diamond)
    for mode in MODES
        xs, ys = series_for(results, mode)
        isempty(ys) && continue
        plot!(plt, xs, ys; label=MODE_LABELS[mode], marker=markers[mode])
    end
    out_path = joinpath(@__DIR__, "leo_thread_scaling_by_mode_$(ROUTE).png")
    savefig(plt, out_path)
    println()
    println("Saved: $(out_path)")
    return plt
end

function main()
    println("Thread-count scaling sweep: route=$(ROUTE), N_SATS=$(N_SATS), modes=$(MODES), thread ladder=$(THREAD_LADDER)")
    println("Estimated wall time: several minutes ($(length(MODES) * length(THREAD_LADDER)) subprocess launches, each with Julia startup + JIT overhead on top of the reported solve time).")
    println()

    results = run_sweep()
    print_summary(results)
    make_plot(results)
end

main()
