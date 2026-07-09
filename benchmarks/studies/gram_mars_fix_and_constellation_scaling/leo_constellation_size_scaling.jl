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
# Julia's thread count is fixed at process startup, and each (N_SATS, mode)
# combination is its own JIT specialization (ODEParams is parameterized on
# both N_sats and the density-model type), so each point is run in a
# separate `julia` subprocess via leo_constellation_size_scaling_worker.jl --
# this avoids accumulating 33 distinct compiled specializations' worth of
# memory in one process (see project findings on this machine's thin memory
# margins). Expect the full sweep to take several minutes.
#
#   julia --project=. benchmarks/studies/leo_constellation_size_scaling.jl
#
# Override the thread count used for every worker subprocess:
#   SPACEAGORA_SCALING_THREADS=8 julia --project=. benchmarks/studies/leo_constellation_size_scaling.jl
#
# Produces two plots in benchmarks/studies/:
#   leo_constellation_size_scaling_with_gram.png
#   leo_constellation_size_scaling_without_gram.png

using Plots

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const WORKER_SCRIPT = joinpath(@__DIR__, "leo_constellation_size_scaling_worker.jl")
const N_SATS_LADDER = [2^k for k in 0:10]  # 1, 2, 4, ..., 1024
const MODES = ("standard", "surrogate", "no_gram")
# Defaults to 4 (not this driver process's own Threads.nthreads(), which is
# silently 1 unless the driver itself is launched with --threads -- the
# worker subprocesses each get their own --threads flag regardless of what
# the driver was started with). 4 matches every other head-to-head timing
# comparison run against this machine (11 physical cores) in this project.
const THREADS = parse(Int, get(ENV, "SPACEAGORA_SCALING_THREADS", "4"))

struct ScalingResult
    n_sats::Int
    mode::String
    median_s::Union{Nothing, Float64}
    output::String
end

function run_worker(n_sats::Int, mode::String)::ScalingResult
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([
        julia_bin, "--project=$(REPO_ROOT)", "--threads=$(THREADS)",
        WORKER_SCRIPT, string(n_sats), mode,
    ])
    io = IOBuffer()
    run(pipeline(cmd; stdout=io, stderr=io); wait=true)
    output = String(take!(io))
    m = match(r"median wall time \(.*?\):\s*([\d.]+)\s*s", output)
    median_s = m === nothing ? nothing : parse(Float64, m.captures[1])
    return ScalingResult(n_sats, mode, median_s, output)
end

function run_sweep()::Vector{ScalingResult}
    results = ScalingResult[]
    total = length(N_SATS_LADDER) * length(MODES)
    done = 0
    for mode in MODES, n_sats in N_SATS_LADDER
        done += 1
        print("[$(done)/$(total)] mode=$(mode) n_sats=$(n_sats) ... ")
        flush(stdout)
        r = run_worker(n_sats, mode)
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

function main()
    println("Constellation-size scaling sweep: N_SATS=$(N_SATS_LADDER), modes=$(MODES), threads=$(THREADS)")
    println("Estimated total wall time: several minutes ($(length(N_SATS_LADDER) * length(MODES)) subprocess launches).")
    println()

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
    make_plots(results)
end

main()
