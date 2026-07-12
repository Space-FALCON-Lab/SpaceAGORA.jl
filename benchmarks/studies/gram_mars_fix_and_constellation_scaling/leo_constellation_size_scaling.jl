# Constellation-size scaling study for a LEO Earth constellation at a fixed
# 600s mission time, comparing:
#
#   "with GRAM"    -- standard (real, no-lookahead-cache) GRAM vs the GRAM
#                     offline surrogate, plotted together
#   "without GRAM" -- NoAtmosphereModel baseline (gravity-only dynamics, no
#                     aerodynamic effector at all), isolating the
#                     N-satellite ODE-integration cost with zero
#                     atmosphere-model overhead
#
# across N_SATS = 1, 2, 4, 8, ..., 1024 (powers of 2).
#
# All 33 (N_SATS, mode) points run in ONE process (leo_constellation_size_scaling_point.jl's
# run_scaling_point), not one subprocess per point. That subprocess-per-point
# design used to be necessary because ODEParams was parameterized on N_sats as
# well as the density-model type, so every distinct satellite count forced a
# fresh JIT specialization of the whole RHS/solver pipeline -- accumulating 33
# specializations' worth of compiled code in one process. N_sats is now a
# runtime field (src/core/types/runtime_types.jl), so every n_sats value for a
# given mode shares one compiled specialization; only the 3 modes still differ
# in type. A single point's exception is caught and recorded as a failure
# (status="failed" in the CSV) rather than aborting the whole sweep, since
# there's no longer a process boundary to contain it.
#
#   julia --project=. --threads=4 benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_constellation_size_scaling.jl
#
# Thread count is fixed at this process's own startup (no more per-worker
# --threads flag to override, since there's no separate worker process) --
# pass --threads=N directly to this invocation instead of via
# SPACEAGORA_SCALING_THREADS.
#
# Produces two plots in benchmarks/studies/:
#   leo_constellation_size_scaling_with_gram.png
#   leo_constellation_size_scaling_without_gram.png
# and the raw per-point results in:
#   leo_constellation_size_scaling_summary.csv

using Plots

include(joinpath(@__DIR__, "leo_constellation_size_scaling_point.jl"))

const N_SATS_LADDER = [2^k for k in 0:10]  # 1, 2, 4, ..., 1024
const MODES = ("standard", "surrogate", "no_gram")

struct ScalingResult
    n_sats::Int
    mode::String
    median_s::Union{Nothing, Float64}
    output::String
end

function run_point(n_sats::Int, mode::String)::ScalingResult
    # invokelatest sidesteps a Julia 1.12 world-age error on first-ever GRAM
    # model construction in this process: GRAMSuite is loaded dynamically
    # (ensure_gramsuite_loaded!'s `@eval import GRAMSuite`) and lazily defines
    # some methods (e.g. `set_library!`) on first use, after this driver's own
    # functions were already compiled against the pre-load world age. Only
    # matters for standard/surrogate (both construct a GRAM model); no_gram
    # never touches GRAMSuite so it was unaffected.
    median_s = try
        Base.invokelatest(run_scaling_point, n_sats, mode)
    catch err
        return ScalingResult(n_sats, mode, nothing, sprint(showerror, err))
    end
    return ScalingResult(n_sats, mode, median_s, "")
end

function run_sweep()::Vector{ScalingResult}
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
            println("$(round(r.median_s; digits=4)) s")
        end
        GC.gc()
    end
    return results
end

function series_for(results::Vector{ScalingResult}, mode::String)
    xs = Int[]
    ys = Float64[]
    for n_sats in N_SATS_LADDER
        r = only(filter(r -> r.mode == mode && r.n_sats == n_sats, results))
        r.median_s === nothing && continue
        push!(xs, n_sats)
        push!(ys, r.median_s)
    end
    return xs, ys
end

# Adds a dashed O(N) reference line anchored at (xs[1], ys[1]) -- on a
# log-log plot, perfect linear scaling from that anchor point is a straight
# line, so deviation of the real curve above/below this line is a visual
# read of super-/sub-linear scaling with constellation size.
function add_linear_reference!(plt, xs, ys)
    isempty(xs) && return plt
    x1, y1 = xs[1], ys[1]
    ref_y = y1 .* (xs ./ x1)
    plot!(plt, xs, ref_y; label="linear (O(N)) reference", linestyle=:dash, color=:gray, marker=:none)
    return plt
end

function make_plots(results::Vector{ScalingResult})
    std_x, std_y = series_for(results, "standard")
    sur_x, sur_y = series_for(results, "surrogate")
    nog_x, nog_y = series_for(results, "no_gram")

    with_gram_plot = plot(
        std_x, std_y;
        label="standard (real GRAM, no lookahead)",
        marker=:circle,
        xscale=:log2,
        yscale=:log10,
        xlabel="Constellation size (satellites)",
        ylabel="Median wall time (s)",
        title="LEO constellation scaling with GRAM\n(Earth, 600s mission)",
        legend=:topleft,
        size=(700, 500),
        margin=5Plots.mm,
    )
    plot!(with_gram_plot, sur_x, sur_y; label="surrogate GRAM", marker=:square)
    add_linear_reference!(with_gram_plot, std_x, std_y)

    without_gram_plot = plot(
        nog_x, nog_y;
        label="no atmosphere model (gravity-only)",
        marker=:circle,
        xscale=:log2,
        yscale=:log10,
        xlabel="Constellation size (satellites)",
        ylabel="Median wall time (s)",
        title="LEO constellation scaling without GRAM\n(Earth, 600s mission)",
        legend=:topleft,
        color=:black,
        size=(700, 500),
        margin=5Plots.mm,
    )
    add_linear_reference!(without_gram_plot, nog_x, nog_y)

    with_path = joinpath(@__DIR__, "leo_constellation_size_scaling_with_gram.png")
    without_path = joinpath(@__DIR__, "leo_constellation_size_scaling_without_gram.png")
    savefig(with_gram_plot, with_path)
    savefig(without_gram_plot, without_path)
    println("Saved: $(with_path)")
    println("Saved: $(without_path)")
    return with_gram_plot, without_gram_plot
end

function save_csv(results::Vector{ScalingResult})
    path = joinpath(@__DIR__, "leo_constellation_size_scaling_summary.csv")
    open(path, "w") do io
        println(io, "n_sats,mode,median_wall_time_s,status")
        for r in results
            status = r.median_s === nothing ? "failed" : "ok"
            value = r.median_s === nothing ? "" : string(r.median_s)
            println(io, "$(r.n_sats),$(r.mode),$(value),$(status)")
        end
    end
    println("Saved: $(path)")
end

function prewarm_gram!()
    # GRAMSuite's package-extension loading can itself trigger a *second*,
    # later world-age bump partway through the very first GRAM/surrogate
    # construction in this process -- after run_point's own invokelatest
    # wrapper has already taken its snapshot, so that one wrapper doesn't
    # cover it (this bit exactly the standard/N=1 point in an earlier run,
    # while every later point was unaffected since the extension only loads
    # once). Force that whole load-and-settle sequence to happen here, before
    # the real (timed, recorded) sweep starts, discarding any result/error.
    for mode in ("standard", "surrogate")
        try
            Base.invokelatest(build_constellation_config, 1, mode)
        catch
        end
    end
    return nothing
end

function main()
    println("Constellation-size scaling sweep: N_SATS=$(N_SATS_LADDER), modes=$(MODES), threads=$(Threads.nthreads())")
    println("Running all $(length(N_SATS_LADDER) * length(MODES)) points in-process (no per-point subprocess).")
    println()

    prewarm_gram!()
    results = run_sweep()

    println()
    println("Summary:")
    println(rpad("n_sats", 10), rpad("standard_s", 14), rpad("surrogate_s", 14), "no_gram_s")
    for n_sats in N_SATS_LADDER
        cells = String[rpad(string(n_sats), 10)]
        for mode in ("standard", "surrogate", "no_gram")
            r = only(filter(r -> r.mode == mode && r.n_sats == n_sats, results))
            push!(cells, rpad(r.median_s === nothing ? "FAILED" : string(round(r.median_s; digits=4)), 14))
        end
        println(join(cells))
    end

    println()
    save_csv(results)
    make_plots(results)
end

main()
