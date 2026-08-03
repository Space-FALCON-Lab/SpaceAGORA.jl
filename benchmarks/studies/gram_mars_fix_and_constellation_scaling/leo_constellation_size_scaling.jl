# Constellation-size x thread-count scaling study for a LEO Earth constellation
# at a fixed 600s mission time, comparing:
#
#   "with GRAM"    -- standard (real, no-lookahead-cache) GRAM vs the GRAM
#                     offline surrogate, plotted together
#   "without GRAM" -- NoAtmosphereModel baseline (gravity-only dynamics, no
#                     aerodynamic effector at all), isolating the
#                     N-satellite ODE-integration cost with zero
#                     atmosphere-model overhead
#
# across N_SATS = 1, 2, 4, 8, ..., 1024 (powers of 2) AND across a thread
# ladder capped at 64 threads (default 1, 2, 4, 8, 16, 32, 64; override with
# SPACEAGORA_SCALING_THREADS as a comma-separated list, e.g. "1,8,64" --
# values above 64 are rejected since this sweep is deliberately scoped to
# "up to 64 threads").
#
#   julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_constellation_size_scaling.jl
#
# Julia's thread count is fixed at process startup, so each thread count in
# the ladder runs as its OWN subprocess (this same script, re-invoked with
# SPACEAGORA_SCALING_THREAD_WORKER=1 and --threads=<N>) -- that's the only
# axis that still needs a process boundary. Within one thread count, all 33
# (n_sats, mode) points run back-to-back in that one subprocess
# (leo_constellation_size_scaling_point.jl's run_scaling_point), not one
# subprocess per point -- ODEParams/SharedBuffers's N_sats is a runtime field
# (src/core/types/runtime_types.jl), so every n_sats value for a given mode
# shares one compiled specialization and only the 3 modes still differ in
# type. A single point's exception is caught and recorded as a failure
# (status="failed" in the CSV) rather than aborting its thread-count worker;
# a whole thread-count worker crashing outright is caught by the top-level
# orchestrator and only drops that thread count's points from the merged
# results, so the rest of the ladder still completes.
#
# Each point also records its own resource footprint -- peak RSS and mean/
# peak CPU%, sampled from THIS process only via `ps -p <pid>` polling
# (resource_monitor.jl) so other load on the machine never pollutes the
# numbers -- across its timed repeats (warmup excluded, since first-use JIT
# cost isn't representative of steady state).
#
# Produces four plots in this directory:
#   leo_constellation_size_scaling_with_gram.png     -- one small-multiple panel
#                                                        per thread count, standard
#                                                        vs surrogate lines
#   leo_constellation_size_scaling_without_gram.png  -- one plot, one line per
#                                                        thread count (light->dark)
#   leo_constellation_size_scaling_resource_ram.png  -- peak RSS, one panel per
#                                                        mode, one line per thread count
#   leo_constellation_size_scaling_resource_cpu.png  -- mean CPU%, same layout
# and the raw per-point results in:
#   leo_constellation_size_scaling_summary.csv

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SELF_PATH = @__FILE__
const N_SATS_LADDER = [2^k for k in 0:10]  # 1, 2, 4, ..., 1024
const MODES = ("standard", "surrogate", "no_gram")
const MODE_LABELS = Dict(
    "standard" => "native GRAM (no lookahead)",
    "surrogate" => "GRAM surrogate",
    "no_gram" => "gravity-only (no atmosphere)",
)
const MODE_COLORS = Dict("standard" => "#4C72B0", "surrogate" => "#DD8452")
const IS_THREAD_WORKER = get(ENV, "SPACEAGORA_SCALING_THREAD_WORKER", "0") == "1"

function parse_thread_ladder(raw::String)::Vector{Int}
    values = sort(unique(parse(Int, strip(tok)) for tok in split(raw, ",")))
    isempty(values) && error("Thread ladder must not be empty.")
    all(v -> v >= 1, values) || error("Thread ladder values must be >= 1, got $(values).")
    maximum(values) <= 64 || error(
        "Thread ladder values must be <= 64 (this sweep is scoped to \"up to 64 threads\"); got max=$(maximum(values))."
    )
    return values
end
const THREAD_LADDER = parse_thread_ladder(get(ENV, "SPACEAGORA_SCALING_THREADS", "1,2,4,8,16,32,64"))

if IS_THREAD_WORKER
    include(joinpath(@__DIR__, "leo_constellation_size_scaling_point.jl"))
else
    using Plots
end

struct ScalingResult
    n_sats::Int
    mode::String
    threads::Int
    median_s::Union{Nothing, Float64}
    peak_rss_mb::Float64
    mean_cpu_pct::Float64
    peak_cpu_pct::Float64
    output::String
end

function run_point(n_sats::Int, mode::String)::ScalingResult
    threads = Threads.nthreads()
    try
        median_s, usage = Base.invokelatest(run_scaling_point, n_sats, mode)
        return ScalingResult(n_sats, mode, threads, median_s, usage.peak_rss_mb, usage.mean_cpu_pct, usage.peak_cpu_pct, "")
    catch err
        return ScalingResult(n_sats, mode, threads, nothing, NaN, NaN, NaN, sprint(showerror, err))
    end
end

function run_worker_sweep()::Vector{ScalingResult}
    results = ScalingResult[]
    total = length(N_SATS_LADDER) * length(MODES)
    done = 0
    for mode in MODES, n_sats in N_SATS_LADDER
        done += 1
        print("[$(done)/$(total)] mode=$(mode) n_sats=$(n_sats) ... ")
        flush(stdout)
        r = run_point(n_sats, mode)
        push!(results, r)
        if r.median_s === nothing
            println("FAILED: $(r.output)")
        else
            println("$(round(r.median_s; digits=4)) s, peak_rss=$(round(r.peak_rss_mb; digits=1)) MB, mean_cpu=$(round(r.mean_cpu_pct; digits=1))%")
        end
        GC.gc()
    end
    return results
end

function worker_csv_path(threads::Int)
    joinpath(@__DIR__, "leo_constellation_size_scaling_summary_threads$(threads).csv")
end

function save_csv(results::Vector{ScalingResult}, path::String)
    open(path, "w") do io
        println(io, "n_sats,mode,threads,median_wall_time_s,peak_rss_mb,mean_cpu_pct,peak_cpu_pct,status")
        for r in results
            status = r.median_s === nothing ? "failed" : "ok"
            median_val = r.median_s === nothing ? "" : string(r.median_s)
            rss_val = isnan(r.peak_rss_mb) ? "" : string(r.peak_rss_mb)
            mean_cpu_val = isnan(r.mean_cpu_pct) ? "" : string(r.mean_cpu_pct)
            peak_cpu_val = isnan(r.peak_cpu_pct) ? "" : string(r.peak_cpu_pct)
            println(io, "$(r.n_sats),$(r.mode),$(r.threads),$(median_val),$(rss_val),$(mean_cpu_val),$(peak_cpu_val),$(status)")
        end
    end
    println("Saved: $(path)")
end

function prewarm_gram!()
    # GRAMSuite's package-extension loading can itself trigger a *second*,
    # later world-age bump partway through the very first GRAM/surrogate
    # construction in this process -- after run_point's own invokelatest
    # wrapper has already taken its snapshot, so that one wrapper doesn't
    # cover it. Force that whole load-and-settle sequence to happen here,
    # before the real (timed, recorded) sweep starts, discarding any
    # result/error.
    for mode in ("standard", "surrogate")
        try
            Base.invokelatest(build_constellation_config, 1, mode)
        catch
        end
    end
    return nothing
end

# Runs as its own subprocess for exactly one thread count: sweeps all
# (n_sats, mode) points in-process at Threads.nthreads() threads and writes
# them to a per-thread-count CSV. No plotting here -- the orchestrator merges
# every thread count's CSV before building the combined plots.
function main_worker()
    threads = Threads.nthreads()
    println("Constellation-size scaling sweep (thread-count worker): N_SATS=$(N_SATS_LADDER), modes=$(MODES), threads=$(threads)")
    println()

    prewarm_gram!()
    results = run_worker_sweep()

    println()
    println("Summary (threads=$(threads)):")
    println(rpad("n_sats", 10), rpad("standard_s", 14), rpad("surrogate_s", 14), "no_gram_s")
    for n_sats in N_SATS_LADDER
        cells = String[rpad(string(n_sats), 10)]
        for mode in MODES
            r = only(filter(r -> r.mode == mode && r.n_sats == n_sats, results))
            push!(cells, rpad(r.median_s === nothing ? "FAILED" : string(round(r.median_s; digits=4)), 14))
        end
        println(join(cells))
    end

    println()
    save_csv(results, worker_csv_path(threads))
end

# Spawns the thread-count worker subprocess for `t` threads, inheriting this
# process's stdout/stderr so the worker's own progress prints stream live.
# Returns true on success; a crashed worker only warns and drops that thread
# count's points from the merged results (mirrors leo_thread_scaling_by_mode.jl's
# handling of a crashed worker) rather than aborting the whole ladder.
function run_thread_count!(t::Int)::Bool
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([julia_bin, "--project=$(REPO_ROOT)", "--threads=$(t)", SELF_PATH])
    println("=== threads=$(t) ===")
    flush(stdout)
    try
        withenv("SPACEAGORA_SCALING_THREAD_WORKER" => "1") do
            run(cmd)
        end
    catch e
        e isa ProcessFailedException || rethrow()
        println("WARNING: threads=$(t) worker subprocess failed (nonzero exit); its points are missing from the merged results.")
        return false
    end
    return true
end

function load_worker_csv(t::Int)::Vector{ScalingResult}
    path = worker_csv_path(t)
    isfile(path) || return ScalingResult[]
    results = ScalingResult[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        isempty(strip(line)) && continue
        cells = split(line, ",")
        n_sats = parse(Int, cells[1])
        mode = cells[2]
        threads = parse(Int, cells[3])
        status = cells[8]
        median_s = status == "ok" ? parse(Float64, cells[4]) : nothing
        peak_rss = isempty(cells[5]) ? NaN : parse(Float64, cells[5])
        mean_cpu = isempty(cells[6]) ? NaN : parse(Float64, cells[6])
        peak_cpu = isempty(cells[7]) ? NaN : parse(Float64, cells[7])
        push!(results, ScalingResult(n_sats, mode, threads, median_s, peak_rss, mean_cpu, peak_cpu, ""))
    end
    rm(path; force=true)
    return results
end

function run_orchestrator_sweep()::Vector{ScalingResult}
    all_results = ScalingResult[]
    for t in THREAD_LADDER
        ok = run_thread_count!(t)
        ok && append!(all_results, load_worker_csv(t))
    end
    return all_results
end

function series_for(results::Vector{ScalingResult}, mode::String, threads::Int)
    xs = Int[]
    ys = Float64[]
    for n_sats in N_SATS_LADDER
        matches = filter(r -> r.mode == mode && r.threads == threads && r.n_sats == n_sats, results)
        isempty(matches) && continue
        r = only(matches)
        r.median_s === nothing && continue
        push!(xs, n_sats)
        push!(ys, r.median_s)
    end
    return xs, ys
end

function resource_series_for(results::Vector{ScalingResult}, mode::String, threads::Int, metric::Symbol)
    xs = Int[]
    ys = Float64[]
    for n_sats in N_SATS_LADDER
        matches = filter(r -> r.mode == mode && r.threads == threads && r.n_sats == n_sats, results)
        isempty(matches) && continue
        r = only(matches)
        val = getfield(r, metric)
        isnan(val) && continue
        push!(xs, n_sats)
        push!(ys, val)
    end
    return xs, ys
end

# Adds a dashed O(N) reference line anchored at (xs[1], ys[1]) -- on a
# log-log plot, perfect linear scaling from that anchor point is a straight
# line, so deviation of the real curve above/below this line is a visual
# read of super-/sub-linear scaling with constellation size.
function add_linear_reference!(plt, xs, ys; label="linear (O(N)) reference")
    isempty(xs) && return plt
    x1, y1 = xs[1], ys[1]
    ref_y = y1 .* (xs ./ x1)
    plot!(plt, xs, ref_y; label=label, linestyle=:dash, color=:gray, marker=:none)
    return plt
end

function facet_layout(n::Int)
    cols = n <= 1 ? 1 : ceil(Int, sqrt(n))
    rows = ceil(Int, n / cols)
    return (rows, cols)
end

# Sequential, single-hue (light -> dark) color ramp keyed by thread count --
# thread count is a magnitude/ordinal dimension here, not a categorical
# identity, so it gets one hue rather than a cycled rainbow. Sampled away
# from the ramp's near-white end so the lowest thread count stays visible
# against a white plot background.
function thread_color_map(threads_list::Vector{Int})
    n = length(threads_list)
    fracs = n == 1 ? [0.85] : collect(range(0.35, 0.95; length=n))
    grad = cgrad(:Blues)
    return Dict(t => grad[f] for (t, f) in zip(threads_list, fracs))
end

function make_plots(results::Vector{ScalingResult})
    rows, cols = facet_layout(length(THREAD_LADDER))

    with_gram_panels = map(THREAD_LADDER) do t
        std_x, std_y = series_for(results, "standard", t)
        sur_x, sur_y = series_for(results, "surrogate", t)
        sp = plot(
            xscale=:log2, yscale=:log10,
            title="$(t) threads", titlefontsize=9,
            xlabel="N sats", ylabel="Median s",
            legend=(t == THREAD_LADDER[1] ? :topleft : false),
            guidefontsize=8, tickfontsize=7, legendfontsize=7,
        )
        isempty(std_y) || plot!(sp, std_x, std_y; label="standard", color=MODE_COLORS["standard"], marker=:circle, markersize=3)
        isempty(sur_y) || plot!(sp, sur_x, sur_y; label="surrogate", color=MODE_COLORS["surrogate"], marker=:square, markersize=3)
        add_linear_reference!(sp, std_x, std_y; label="")
        sp
    end
    with_gram_plot = plot(
        with_gram_panels...; layout=(rows, cols), size=(320 * cols, 260 * rows + 40),
        plot_title="LEO constellation scaling with GRAM, by thread count (Earth, 600s mission)", plot_titlefontsize=11,
        margin=3Plots.mm,
    )

    colors = thread_color_map(THREAD_LADDER)
    without_gram_plot = plot(
        xlabel="Constellation size (satellites)", ylabel="Median wall time (s)",
        title="LEO constellation scaling without GRAM, by thread count\n(Earth, 600s mission)",
        xscale=:log2, yscale=:log10, legend=:topleft, size=(750, 550), margin=5Plots.mm,
    )
    for t in THREAD_LADDER
        xs, ys = series_for(results, "no_gram", t)
        isempty(ys) && continue
        plot!(without_gram_plot, xs, ys; label="$(t) threads", color=colors[t], marker=:circle, markersize=4)
    end
    ref_xs, ref_ys = series_for(results, "no_gram", THREAD_LADDER[1])
    add_linear_reference!(without_gram_plot, ref_xs, ref_ys)

    with_path = joinpath(@__DIR__, "leo_constellation_size_scaling_with_gram.png")
    without_path = joinpath(@__DIR__, "leo_constellation_size_scaling_without_gram.png")
    savefig(with_gram_plot, with_path)
    savefig(without_gram_plot, without_path)
    println("Saved: $(with_path)")
    println("Saved: $(without_path)")
    return with_gram_plot, without_gram_plot
end

function make_resource_plot(results::Vector{ScalingResult}, metric::Symbol, ylabel::String, title::String, out_name::String; log_y::Bool=true)
    colors = thread_color_map(THREAD_LADDER)
    panels = map(MODES) do mode
        sp = plot(
            xscale=:log2, yscale=(log_y ? :log10 : :identity),
            title=MODE_LABELS[mode], titlefontsize=9,
            xlabel="N sats", ylabel=ylabel,
            legend=(mode == MODES[1] ? :topleft : false),
            guidefontsize=8, tickfontsize=7, legendfontsize=7,
        )
        for t in THREAD_LADDER
            xs, ys = resource_series_for(results, mode, t, metric)
            isempty(ys) && continue
            plot!(sp, xs, ys; label="$(t) threads", color=colors[t], marker=:circle, markersize=3)
        end
        sp
    end
    plt = plot(
        panels...; layout=(1, length(MODES)), size=(350 * length(MODES), 420),
        plot_title=title, plot_titlefontsize=11, margin=3Plots.mm,
    )
    out_path = joinpath(@__DIR__, out_name)
    savefig(plt, out_path)
    println("Saved: $(out_path)")
    return plt
end

function make_resource_plots(results::Vector{ScalingResult})
    make_resource_plot(
        results, :peak_rss_mb, "Peak RSS (MB)",
        "LEO constellation resource footprint: peak RAM, by mode and thread count (Earth, 600s mission)",
        "leo_constellation_size_scaling_resource_ram.png"; log_y=true,
    )
    make_resource_plot(
        results, :mean_cpu_pct, "Mean CPU (%, this process only)",
        "LEO constellation resource footprint: mean CPU usage, by mode and thread count (Earth, 600s mission)",
        "leo_constellation_size_scaling_resource_cpu.png"; log_y=false,
    )
end

function main_orchestrator()
    n_points = length(N_SATS_LADDER) * length(MODES)
    println("Constellation-size x thread-count scaling sweep: N_SATS=$(N_SATS_LADDER), modes=$(MODES), threads=$(THREAD_LADDER)")
    println("Each of the $(length(THREAD_LADDER)) thread counts runs as its own subprocess of $(n_points) in-process (n_sats, mode) points.")
    println()

    results = run_orchestrator_sweep()

    println()
    println("Summary:")
    for t in THREAD_LADDER
        println()
        println("-- threads=$(t) --")
        println(rpad("n_sats", 10), rpad("standard_s", 14), rpad("surrogate_s", 14), "no_gram_s")
        for n_sats in N_SATS_LADDER
            cells = String[rpad(string(n_sats), 10)]
            for mode in MODES
                matches = filter(r -> r.mode == mode && r.threads == t && r.n_sats == n_sats, results)
                cell = isempty(matches) ? "--" : (only(matches).median_s === nothing ? "FAILED" : string(round(only(matches).median_s; digits=4)))
                push!(cells, rpad(cell, 14))
            end
            println(join(cells))
        end
    end

    println()
    save_csv(results, joinpath(@__DIR__, "leo_constellation_size_scaling_summary.csv"))
    make_plots(results)
    make_resource_plots(results)
end

if IS_THREAD_WORKER
    main_worker()
else
    main_orchestrator()
end
