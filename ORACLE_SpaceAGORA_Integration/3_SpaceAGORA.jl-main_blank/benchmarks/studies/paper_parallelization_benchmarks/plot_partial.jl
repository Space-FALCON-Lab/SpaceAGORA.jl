# Plots the results of a paper-benchmarks run that may still be in progress.
#
# ppc_run_controller only writes its combined raw/aggregated CSVs once every
# worker subprocess for a phase (or, for B4's worker-ladder sweep, a single
# workers_NN sub-run) has finished. Each worker still writes its own row to
# <phase_dir>/worker_rows/perf_*.csv the moment it completes, so this script
# reassembles a snapshot from whatever is on disk right now: finished
# per-phase raw CSVs where available, and raw worker_rows/*.csv files for
# any phase/sub-run that hasn't finished yet.
#
# Usage:
#   julia --project=. benchmarks/studies/paper_parallelization_benchmarks/plot_partial.jl [run_dir]
#
# run_dir defaults to the most recently modified run under
# SPACEAGORA_PPB_OUTDIR / the default paper-benchmarks output root. Output
# (raw CSV, aggregated CSV, plots, short report) is written to
# <run_dir>/partial_snapshot/, overwritten on every invocation so it can be
# re-run repeatedly while the benchmark is still executing.

include(joinpath(@__DIR__, "..", "parallelization_performance", "cli.jl"))
include(joinpath(@__DIR__, "..", "parallelization_performance", "reporting.jl"))
include(joinpath(@__DIR__, "cli.jl"))
include(joinpath(@__DIR__, "reporting.jl"))

function _ppb_latest_run_dir(root::String)::String
    isdir(root) || throw(ArgumentError("No paper-benchmarks output directory found at $(root)."))
    candidates = [
        joinpath(root, d) for d in readdir(root)
        if isdir(joinpath(root, d)) && !startswith(d, "dry_run_")
    ]
    isempty(candidates) && throw(ArgumentError("No runs found under $(root)."))
    return argmax(mtime, candidates)
end

function _ppb_write_partial_report(
    path::String, run_dir::String, collected::NamedTuple, agg::DataFrame
)::Nothing
    open(path, "w") do io
        println(io, "# Partial Snapshot — Paper Parallelization Benchmarks")
        println(io)
        println(io, "Generated: $(now(UTC)) (snapshot of an in-progress or completed run)")
        println(io, "Run directory: $(run_dir)")
        println(io)
        println(io, "## Complete phase/sub-run directories ($(length(collected.complete_dirs)))")
        for d in sort(collected.complete_dirs)
            println(io, "- $(relpath(d, run_dir))")
        end
        println(io)
        println(io, "## In-progress phase/sub-run directories ($(length(collected.partial_dirs)))")
        if isempty(collected.partial_dirs)
            println(io, "_None — every directory found so far has a completed raw CSV._")
        else
            for d in sort(collected.partial_dirs)
                println(io, "- $(relpath(d, run_dir))")
            end
        end
        println(io)
        println(io, "Rows collected: $(nrow(collected.raw)) raw, $(nrow(agg)) aggregated groups.")
    end
    return nothing
end

function main_plot_partial(args::Vector{String}=ARGS)
    run_dir = isempty(args) ? _ppb_latest_run_dir(PPB_DEFAULT_OUTDIR) : abspath(args[1])
    isdir(run_dir) || throw(ArgumentError("Run directory does not exist: $(run_dir)"))
    println("[plot-partial] run_dir = $(run_dir)")

    collected = _ppb_collect_partial(run_dir)
    println("[plot-partial] complete phase/sub-run dirs   = $(length(collected.complete_dirs))")
    println("[plot-partial] in-progress phase/sub-run dirs = $(length(collected.partial_dirs))")
    for d in collected.partial_dirs
        println("  in progress: $(relpath(d, run_dir))")
    end

    if nrow(collected.raw) == 0
        println("[plot-partial] No data found yet — no worker CSVs have been written.")
        return
    end

    agg = _ppb_aggregate(collected.raw)
    if nrow(agg) == 0
        println("[plot-partial] Rows found but none successful yet; nothing to plot.")
        return
    end

    snapshot_dir = joinpath(run_dir, "partial_snapshot")
    mkpath(snapshot_dir)

    raw_path = joinpath(snapshot_dir, "partial_raw.csv")
    agg_path = joinpath(snapshot_dir, "partial_aggregated.csv")
    CSV.write(raw_path, collected.raw)
    CSV.write(agg_path, agg)
    println("[plot-partial] raw CSV        = $(raw_path)")
    println("[plot-partial] aggregated CSV = $(agg_path)")

    plot_paths = _ppb_write_plots(snapshot_dir, agg)
    println("[plot-partial] plots (n=$(length(plot_paths))):")
    for p in plot_paths
        println("  $(p)")
    end

    report_path = joinpath(snapshot_dir, "partial_report.md")
    _ppb_write_partial_report(report_path, run_dir, collected, agg)
    println("[plot-partial] report         = $(report_path)")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_plot_partial()
end
