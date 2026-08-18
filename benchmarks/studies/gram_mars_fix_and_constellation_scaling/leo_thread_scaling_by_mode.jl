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
# Each point also records the worker subprocess's own resource footprint --
# peak RSS and mean/peak CPU%, sampled by pid via `ps` (resource_monitor.jl)
# so other load on the machine never pollutes the numbers -- plotted
# separately from the wall-time curve.
#
#   julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_thread_scaling_by_mode.jl
#
# Default thread ladder now goes up to 64 threads (1, 2, 4, 8, 16, 32, 64).
# Override satellite count / thread ladder:
#   SPACEAGORA_TS_N_SATS=512 SPACEAGORA_TS_THREADS=1,2,4,8,16,32,64 julia --project=. \
#       benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_thread_scaling_by_mode.jl

using Plots

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WORKER_SCRIPT = joinpath(@__DIR__, "leo_constellation_size_scaling_worker.jl")
include(joinpath(@__DIR__, "resource_monitor.jl"))
const N_SATS = parse(Int, get(ENV, "SPACEAGORA_TS_N_SATS", "512"))
const THREAD_LADDER = [parse(Int, strip(tok)) for tok in split(get(ENV, "SPACEAGORA_TS_THREADS", "1,2,4,8,16,32,64"), ",")]
maximum(THREAD_LADDER) <= 64 || error(
    "Thread ladder values must be <= 64 (this sweep is scoped to \"up to 64 threads\"); got max=$(maximum(THREAD_LADDER))."
)
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
const MARKERS = Dict("no_gram" => :circle, "surrogate" => :square, "surrogate_on" => :utriangle, "standard" => :diamond)

struct ThreadResult
    mode::String
    threads::Int
    median_s::Union{Nothing, Float64}
    peak_rss_mb::Float64
    mean_cpu_pct::Float64
    peak_cpu_pct::Float64
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
    # concurrency *within* one coupled solve). wait=false (rather than the previous
    # wait=true + catch ProcessFailedException) so the worker's own pid is available
    # to resource_monitor.jl for the whole subprocess's lifetime; a crashed/nonzero-exit
    # worker is detected via `success(proc)` afterwards instead of an exception, and
    # only drops that point (median_s stays nothing) rather than aborting the whole
    # sweep so the rest of the ladder/modes still complete.
    usage = _EMPTY_RESOURCE_USAGE
    withenv(env_overrides...) do
        proc = run(pipeline(cmd; stdout=io, stderr=io); wait=false)
        mon = start_resource_monitor(getpid(proc))
        wait(proc)
        usage = stop_and_collect!(mon)
        success(proc) || println("  (worker exited nonzero for mode=$(mode) threads=$(threads); see raw output on failure below)")
    end
    output = String(take!(io))
    m = match(r"median wall time \(.*?\):\s*([\d.]+)\s*s", output)
    median_s = m === nothing ? nothing : parse(Float64, m.captures[1])
    return ThreadResult(mode, threads, median_s, usage.peak_rss_mb, usage.mean_cpu_pct, usage.peak_cpu_pct, output)
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
            println("$(round(r.median_s; digits=4)) s, peak_rss=$(round(r.peak_rss_mb; digits=1)) MB, mean_cpu=$(round(r.mean_cpu_pct; digits=1))%")
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

function resource_series_for(results::Vector{ThreadResult}, mode::String, metric::Symbol)
    xs = Int[]
    ys = Float64[]
    for threads in THREAD_LADDER
        r = only(filter(r -> r.mode == mode && r.threads == threads, results))
        val = getfield(r, metric)
        isnan(val) && continue
        push!(xs, threads)
        push!(ys, val)
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
        println(rpad("threads", 10), rpad("median_s", 12), rpad("speedup", 10), rpad("efficiency", 12), rpad("rss_MB", 10), "cpu_%")
        for t in xs
            r = only(filter(r -> r.mode == mode && r.threads == t, results))
            y = r.median_s
            speedup = baseline / y
            efficiency = speedup / t
            rss_str = isnan(r.peak_rss_mb) ? "--" : string(round(r.peak_rss_mb; digits=1))
            cpu_str = isnan(r.mean_cpu_pct) ? "--" : string(round(r.mean_cpu_pct; digits=1))
            println(
                rpad(string(t), 10),
                rpad(string(round(y; digits=4)), 12),
                rpad(string(round(speedup; digits=3)), 10),
                rpad(string(round(efficiency; digits=3)), 12),
                rpad(rss_str, 10),
                cpu_str,
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
    for mode in MODES
        xs, ys = series_for(results, mode)
        isempty(ys) && continue
        plot!(plt, xs, ys; label=MODE_LABELS[mode], marker=MARKERS[mode])
    end
    out_path = joinpath(@__DIR__, "leo_thread_scaling_by_mode_$(ROUTE).png")
    savefig(plt, out_path)
    println()
    println("Saved: $(out_path)")
    return plt
end

# Separate figure (not overlaid on the wall-time plot, per "one axis" -- RSS
# and CPU% are a different unit/scale than wall time) with peak RSS and mean
# CPU% side by side, one line per mode. Both metrics come from resource_monitor.jl
# sampling the worker subprocess by pid, so they reflect only that process --
# other load on the machine doesn't leak into either panel.
function make_resource_plot(results::Vector{ThreadResult})
    rss_panel = plot(
        xlabel="Threads", ylabel="Peak RSS (MB)", title="Peak RAM",
        xscale=:log2, legend=:topleft, guidefontsize=9, titlefontsize=10,
    )
    cpu_panel = plot(
        xlabel="Threads", ylabel="Mean CPU (%, this process only)", title="Mean CPU usage",
        xscale=:log2, legend=false, guidefontsize=9, titlefontsize=10,
    )
    for mode in MODES
        rss_xs, rss_ys = resource_series_for(results, mode, :peak_rss_mb)
        cpu_xs, cpu_ys = resource_series_for(results, mode, :mean_cpu_pct)
        isempty(rss_ys) || plot!(rss_panel, rss_xs, rss_ys; label=MODE_LABELS[mode], marker=MARKERS[mode])
        isempty(cpu_ys) || plot!(cpu_panel, cpu_xs, cpu_ys; label=MODE_LABELS[mode], marker=MARKERS[mode])
    end
    plt = plot(
        rss_panel, cpu_panel; layout=(1, 2), size=(1100, 500),
        plot_title="Resource footprint by density-model mode ($(ROUTE))\n(LEO, $(N_SATS) sats, 600s mission)",
        plot_titlefontsize=11, margin=6Plots.mm,
    )
    out_path = joinpath(@__DIR__, "leo_thread_scaling_by_mode_resource_$(ROUTE).png")
    savefig(plt, out_path)
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
    make_resource_plot(results)
end

main()
