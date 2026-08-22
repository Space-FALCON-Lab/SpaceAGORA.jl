using Statistics

# ── Helpers ──────────────────────────────────────────────────────────────────

function _ppb_n_sat(case_name::String)::Int
    m = match(r"(\d+)sat", case_name)
    m !== nothing && return parse(Int, m.captures[1])
    occursin("single", case_name) && return 1
    m = match(r"(?:multi|callback)_(\d+)_", case_name)
    m !== nothing && return parse(Int, m.captures[1])
    return 1
end

function _ppb_phase_from_path(dirpath::String, root::String)::String
    rel = relpath(dirpath, root)
    rel == "." && return "unknown"
    return split(rel, Base.Filesystem.path_separator)[1]
end

# ── Data collection ───────────────────────────────────────────────────────────

# `ppc_run_controller` only writes its combined `parallelization_performance_raw_*.csv`
# once every worker subprocess for that controller call (a phase, or for B4/B6-style
# worker-ladder sweeps, a single workers_NN sub-run) has finished. Each worker still
# writes its own row to <phase_dir>/worker_rows/perf_*.csv the moment it completes, so
# collection here reassembles a snapshot from whatever is on disk right now: finished
# per-phase raw CSVs where available, and raw worker_rows/*.csv files for any
# phase/sub-run that hasn't finished yet. This makes it safe to call repeatedly while
# the benchmark is still executing, not just once at the very end.

function _ppb_read_csv_safe(path::String)::DataFrame
    try
        return ppc_read_optional(path)
    catch err
        @warn "Could not read $(path) (likely still being written); skipping for now" exception=err
        return DataFrame()
    end
end

function _ppb_partial_worker_csvs(worker_rows_dir::String)::DataFrame
    frames = DataFrame[]
    for file in readdir(worker_rows_dir)
        (startswith(file, "perf_") && endswith(file, ".csv")) || continue
        df = _ppb_read_csv_safe(joinpath(worker_rows_dir, file))
        nrow(df) > 0 && push!(frames, df)
    end
    isempty(frames) && return DataFrame()
    return vcat(frames...; cols=:union)
end

"""
    _ppb_collect_partial(root_dir) -> (raw, complete_dirs, partial_dirs)

Walk `root_dir` and assemble a raw results DataFrame from whatever is on
disk: a finished `parallelization_performance_raw_*.csv` where a phase (or
worker-ladder sub-run) has completed, otherwise whatever individual
`worker_rows/perf_*.csv` files that sub-run has produced so far.
"""
function _ppb_collect_partial(root_dir::String)::NamedTuple
    isdir(root_dir) || return (raw=DataFrame(), complete_dirs=String[], partial_dirs=String[])

    frames = DataFrame[]
    complete_dirs = String[]
    partial_dirs = String[]
    claimed = Set{String}()

    for (dirpath, _, files) in walkdir(root_dir)
        finished = [f for f in files if startswith(f, "parallelization_performance_raw_") && endswith(f, ".csv")]
        isempty(finished) && continue
        for file in finished
            df = _ppb_read_csv_safe(joinpath(dirpath, file))
            nrow(df) == 0 && continue
            df[!, :phase_id] .= _ppb_phase_from_path(dirpath, root_dir)
            df[!, :run_status] .= "complete"
            push!(frames, df)
        end
        push!(complete_dirs, dirpath)
        push!(claimed, dirpath)
    end

    for (dirpath, dirs, _) in walkdir(root_dir)
        dirpath in claimed && continue
        "worker_rows" in dirs || continue
        df = _ppb_partial_worker_csvs(joinpath(dirpath, "worker_rows"))
        nrow(df) == 0 && continue
        df[!, :phase_id] .= _ppb_phase_from_path(dirpath, root_dir)
        df[!, :run_status] .= "partial"
        push!(frames, df)
        push!(partial_dirs, dirpath)
    end

    raw = isempty(frames) ? DataFrame() : vcat(frames...; cols=:union)
    # CSV.read infers compact InlineStrings types (e.g. String31), but
    # _ppb_n_sat and friends are typed on plain String; normalize here so
    # _ppb_aggregate doesn't choke on a partial snapshot's column types.
    "case" in names(raw) && (raw[!, :case] = String.(raw.case))
    return (raw=raw, complete_dirs=complete_dirs, partial_dirs=partial_dirs)
end

# ── Router regret (point 8: report router regret relative to the best tested
# static route, (T_selected - T_best) / T_best) ────────────────────────────────

# "Static" routes are fixed, non-adaptive parallel profiles (R0-R3); "adaptive"
# routes let the runtime pick the schedule per call/workload (R4-R5). Regret is
# computed for the adaptive routes against the best static route measured on
# the *same* workload point (phase_id, case, thread_count, process_workers,
# mc_samples) -- i.e. what the adaptive router actually cost you relative to
# the best fixed route you could have picked for that exact point, matching
# the review's point 8 formula.
const PPB_STATIC_MODES   = Set(["serial", "outer_threads", "outer_process", "inner_only", "outer_inner_static"])
const PPB_ADAPTIVE_MODES = Set(["outer_inner_adaptive", "full_smart"])

function _ppb_add_router_regret!(agg::DataFrame)
    agg[!, :regret_vs_best_static] = Vector{Union{Missing, Float64}}(missing, nrow(agg))
    nrow(agg) == 0 && return agg
    point_key = [k for k in [:phase_id, :case, :thread_count, :process_workers, :mc_samples] if string(k) in names(agg)]
    isempty(point_key) && return agg

    for sub in groupby(agg, point_key)
        static_times = [t for (m, t) in zip(sub.mode, sub.wall_time_median_s) if m in PPB_STATIC_MODES && !ismissing(t) && t > 0.0]
        isempty(static_times) && continue
        best_static = minimum(static_times)
        for i in 1:nrow(sub)
            t = sub.wall_time_median_s[i]
            (ismissing(t) || t <= 0.0) && continue
            sub.regret_vs_best_static[i] = (t - best_static) / best_static
        end
    end
    return agg
end

# ── Aggregation ───────────────────────────────────────────────────────────────

function _ppb_aggregate(raw::DataFrame)::DataFrame
    nrow(raw) == 0 && return DataFrame()
    cols = names(raw)

    group_keys = [:phase_id, :case, :mode, :thread_count, :process_workers, :mc_samples]
    group_keys = [k for k in group_keys if string(k) in cols]

    # Restrict to successful rows.
    ok = "success" in cols ? coalesce.(raw.success, false) : fill(true, nrow(raw))
    df = raw[ok, :]
    nrow(df) == 0 && return DataFrame()

    agg = combine(
        groupby(df, group_keys),
        :wall_time_s => (x -> median(collect(skipmissing(x)))) => :wall_time_median_s,
        :wall_time_s => (x -> mean(collect(skipmissing(x))))   => :wall_time_mean_s,
        :wall_time_s => (x -> std(collect(skipmissing(x))))    => :wall_time_std_s,
        :wall_time_s => length                                   => :n_repeats,
        :throughput_samples_per_s => (x -> mean(collect(skipmissing(x)))) => :throughput_mean,
    )

    # Serial baseline: join within (phase_id, case, mc_samples, process_workers) so
    # each parallel row is compared against the serial run from the same sub-run.
    serial_key = [k for k in [:phase_id, :case, :mc_samples, :process_workers] if k in Symbol.(names(agg))]
    serial_df  = agg[agg.mode .== "serial", [serial_key..., :wall_time_median_s]]
    rename!(serial_df, :wall_time_median_s => :serial_median_s)
    agg = leftjoin(agg, serial_df; on=serial_key)

    agg[!, :speedup] = [
        (ismissing(s) || t <= 0.0) ? missing : s / t
        for (s, t) in zip(agg.serial_median_s, agg.wall_time_median_s)
    ]
    agg[!, :thread_efficiency] = [
        (ismissing(sp) || tc <= 1) ? missing : sp / max(1, tc)
        for (sp, tc) in zip(agg.speedup, agg.thread_count)
    ]
    if "process_workers" in names(agg)
        agg[!, :process_efficiency] = [
            (ismissing(sp) || pw <= 1) ? missing : sp / max(1, pw)
            for (sp, pw) in zip(agg.speedup, agg.process_workers)
        ]
    end
    agg[!, :n_sat] = _ppb_n_sat.(agg.case)
    _ppb_add_router_regret!(agg)
    return agg
end

# ── Individual plots ──────────────────────────────────────────────────────────

function _ppb_plot_constellation_scaling(
    agg::DataFrame, phase_id::String, title::String, filename::String, outdir::String
)::String
    df = agg[coalesce.(agg.phase_id, "") .== phase_id, :]
    nrow(df) == 0 && return ""
    df = df[df.mode .!= "serial", :]
    nrow(df) == 0 && return ""

    max_t   = maximum(df.thread_count)
    plot_df = df[df.thread_count .== max_t, :]
    nrow(plot_df) == 0 && return ""

    modes     = sort(unique(plot_df.mode))
    n_sat_ref = sort(unique(plot_df.n_sat))

    p = Plots.plot(
        title       = title,
        xlabel      = "Spacecraft (N)",
        ylabel      = "Speedup vs. serial (×)",
        xscale      = :log10,
        legend      = :topleft,
        tickfont    = Plots.font(9),
        guidefont   = Plots.font(11),
        titlefont   = Plots.font(11),
        legendfont  = Plots.font(9),
        size        = (700, 480),
        left_margin = 12Plots.PlotMeasures.mm,
        bottom_margin = 8Plots.PlotMeasures.mm,
    )

    # Ideal linear reference capped at max_t.
    n_ref = sort(unique(vcat(n_sat_ref, [1, max_t])))
    Plots.plot!(p, n_ref, min.(n_ref, max_t);
        label="ideal (capped at $(max_t) threads)",
        linestyle=:dash, color=:grey, linewidth=1)

    for mode in modes
        sub = plot_df[plot_df.mode .== mode, :]
        nrow(sub) == 0 && continue
        sub = sort(sub, :n_sat)
        Plots.plot!(p, sub.n_sat, coalesce.(sub.speedup, NaN);
            label=mode, marker=:circle, linewidth=2)
    end

    path = joinpath(outdir, filename)
    Plots.savefig(p, path)
    return path
end

function _ppb_plot_thread_scaling(agg::DataFrame, outdir::String)::String
    # Speedup vs thread count for B1 and B2, outer_threads mode, large N_sat values.
    phases   = ["B1", "B2"]
    target_n = [64, 256, 1024]
    mode_sel = "outer_threads"

    df = agg[
        (coalesce.(agg.phase_id, "") .∈ Ref(phases)) .&
        (agg.mode .== mode_sel) .&
        (agg.n_sat .∈ Ref(target_n)),
        :
    ]
    nrow(df) == 0 && return ""

    p = Plots.plot(
        title       = "Thread Scaling — outer_threads",
        xlabel      = "Thread count",
        ylabel      = "Speedup vs. serial (×)",
        legend      = :topleft,
        tickfont    = Plots.font(9),
        guidefont   = Plots.font(11),
        titlefont   = Plots.font(11),
        legendfont  = Plots.font(9),
        size        = (700, 480),
        left_margin = 12Plots.PlotMeasures.mm,
        bottom_margin = 8Plots.PlotMeasures.mm,
    )

    for phase in phases
        for n in sort(target_n)
            sub = df[(coalesce.(df.phase_id, "") .== phase) .& (df.n_sat .== n), :]
            nrow(sub) == 0 && continue
            sub = sort(sub, :thread_count)
            label = "$(phase) / N=$(n)"
            Plots.plot!(p, sub.thread_count, coalesce.(sub.speedup, NaN);
                label=label, marker=:circle, linewidth=2)
        end
    end

    # Ideal reference.
    max_t = maximum(df.thread_count)
    Plots.plot!(p, 1:max_t, 1:max_t;
        label="ideal", linestyle=:dash, color=:grey, linewidth=1)

    path = joinpath(outdir, "thread_scaling_outer_threads.png")
    Plots.savefig(p, path)
    return path
end

function _ppb_plot_atmosphere_runtime(agg::DataFrame, outdir::String)::String
    df = agg[coalesce.(agg.phase_id, "") .== "B3", :]
    nrow(df) == 0 && return ""
    df = sort(df, [:case, :mode])

    cases = sort(unique(df.case))
    modes = sort(unique(df.mode))
    n_c, n_m = length(cases), length(modes)
    n_c == 0 || n_m == 0 && return ""

    # Build matrix: rows = cases, columns = modes.
    y = [begin
        sub = df[(df.case .== c) .& (df.mode .== m), :]
        nrow(sub) > 0 ? coalesce(sub.wall_time_median_s[1], NaN) : NaN
    end for c in cases, m in modes]

    p = Plots.bar(
        cases, y;
        label       = permutedims(modes),
        ylabel      = "Median wall time (s)",
        title       = "B3 — Atmosphere Runtime by Mode",
        xrotation   = 20,
        tickfont    = Plots.font(8),
        guidefont   = Plots.font(11),
        titlefont   = Plots.font(11),
        legendfont  = Plots.font(9),
        legend      = :topright,
        size        = (700, 480),
        left_margin = 14Plots.PlotMeasures.mm,
        bottom_margin = 24Plots.PlotMeasures.mm,
    )

    path = joinpath(outdir, "atmosphere_runtime_comparison.png")
    Plots.savefig(p, path)
    return path
end

function _ppb_plot_mc_throughput(agg::DataFrame, outdir::String)::String
    df = agg[(coalesce.(agg.phase_id, "") .== "B4") .& (agg.mode .== "outer_process"), :]
    nrow(df) == 0 && return ""

    p = Plots.plot(
        title       = "B4 — Monte Carlo Throughput vs. Process Workers",
        xlabel      = "Process workers",
        ylabel      = "Speedup vs. serial (×)",
        legend      = :topleft,
        tickfont    = Plots.font(9),
        guidefont   = Plots.font(11),
        titlefont   = Plots.font(11),
        legendfont  = Plots.font(9),
        size        = (700, 480),
        left_margin = 12Plots.PlotMeasures.mm,
        bottom_margin = 8Plots.PlotMeasures.mm,
    )

    # Ideal linear reference.
    max_w = maximum(skipmissing(df.process_workers))
    Plots.plot!(p, 1:max_w, 1:max_w;
        label="ideal", linestyle=:dash, color=:grey, linewidth=1)

    for case in sort(unique(df.case))
        for mc in sort(unique(df.mc_samples))
            sub = df[(df.case .== case) .& (df.mc_samples .== mc), :]
            nrow(sub) == 0 && continue
            sub = sort(sub, :process_workers)
            Plots.plot!(p, sub.process_workers, coalesce.(sub.speedup, NaN);
                label="$(case) / mc=$(mc)", marker=:circle, linewidth=2)
        end
    end

    path = joinpath(outdir, "mc_throughput_scaling.png")
    Plots.savefig(p, path)
    return path
end

# B7's counterpart to _ppb_plot_thread_scaling. Kept separate rather than folded
# into that function's phase list because B7 varies the *case* (four different
# workload shapes) rather than N_sat within one physics configuration, so the
# series have to be keyed on case name -- and because the point of the plot is
# the contrast between where each case's curve flattens, which needs them on one
# set of axes against the ideal line.
function _ppb_plot_heavy_thread_scaling(agg::DataFrame, outdir::String)::String
    df = agg[(coalesce.(agg.phase_id, "") .== "B7") .& (agg.mode .== "outer_threads"), :]
    nrow(df) == 0 && return ""

    p = Plots.plot(
        title       = "B7 — Heavy Constellation Thread Scaling (outer_threads)",
        xlabel      = "Thread count (physical cores)",
        ylabel      = "Speedup vs. serial (×)",
        legend      = :topleft,
        tickfont    = Plots.font(9),
        guidefont   = Plots.font(11),
        titlefont   = Plots.font(11),
        legendfont  = Plots.font(8),
        size        = (700, 480),
        left_margin = 12Plots.PlotMeasures.mm,
        bottom_margin = 8Plots.PlotMeasures.mm,
    )

    for case in sort(unique(df.case))
        sub = sort(df[df.case .== case, :], :thread_count)
        nrow(sub) == 0 && continue
        Plots.plot!(p, sub.thread_count, coalesce.(sub.speedup, NaN);
            label=case, marker=:circle, linewidth=2)
    end

    max_t = maximum(df.thread_count)
    Plots.plot!(p, 1:max_t, 1:max_t;
        label="ideal", linestyle=:dash, color=:grey, linewidth=1)

    path = joinpath(outdir, "b7_heavy_thread_scaling.png")
    Plots.savefig(p, path)
    return path
end

function _ppb_plot_heavy_mc_throughput(agg::DataFrame, outdir::String)::String
    df = agg[(coalesce.(agg.phase_id, "") .== "B8") .& (agg.mode .== "outer_process"), :]
    nrow(df) == 0 && return ""

    p = Plots.plot(
        title       = "B8 — Heavy Monte Carlo Throughput vs. Process Workers",
        xlabel      = "Process workers",
        ylabel      = "Speedup vs. serial (×)",
        legend      = :topleft,
        tickfont    = Plots.font(9),
        guidefont   = Plots.font(11),
        titlefont   = Plots.font(11),
        legendfont  = Plots.font(9),
        size        = (700, 480),
        left_margin = 12Plots.PlotMeasures.mm,
        bottom_margin = 8Plots.PlotMeasures.mm,
    )

    max_w = maximum(skipmissing(df.process_workers))
    Plots.plot!(p, 1:max_w, 1:max_w;
        label="ideal", linestyle=:dash, color=:grey, linewidth=1)

    for mc in sort(unique(df.mc_samples))
        sub = sort(df[df.mc_samples .== mc, :], :process_workers)
        nrow(sub) == 0 && continue
        Plots.plot!(p, sub.process_workers, coalesce.(sub.speedup, NaN);
            label="mc=$(mc)", marker=:circle, linewidth=2)
    end

    path = joinpath(outdir, "b8_heavy_mc_throughput.png")
    Plots.savefig(p, path)
    return path
end

function _ppb_plot_profile_comparison(agg::DataFrame, outdir::String)::String
    df = agg[coalesce.(agg.phase_id, "") .== "B5", :]
    nrow(df) == 0 && return ""

    cases      = sort(unique(df.case))
    mode_order = ["serial", "outer_threads", "inner_only", "outer_inner_static", "outer_inner_adaptive", "full_smart"]
    modes      = [m for m in mode_order if m in unique(df.mode)]
    n_c, n_m   = length(cases), length(modes)
    n_c == 0 || n_m == 0 && return ""

    # Speedup matrix: rows = cases, columns = modes.
    y = [begin
        sub = df[(df.case .== c) .& (df.mode .== m), :]
        nrow(sub) > 0 ? coalesce(sub.speedup[1], NaN) : NaN
    end for c in cases, m in modes]

    p = Plots.bar(
        cases, y;
        label         = permutedims(modes),
        ylabel        = "Speedup vs. serial (×)",
        title         = "B5 — Routing Profile Comparison",
        xrotation     = 20,
        tickfont      = Plots.font(8),
        guidefont     = Plots.font(11),
        titlefont     = Plots.font(11),
        legendfont    = Plots.font(9),
        legend        = :topright,
        size          = (800, 500),
        left_margin   = 14Plots.PlotMeasures.mm,
        bottom_margin = 24Plots.PlotMeasures.mm,
        right_margin  = 32Plots.PlotMeasures.mm,
    )

    path = joinpath(outdir, "profile_comparison.png")
    Plots.savefig(p, path)
    return path
end

function _ppb_plot_router_regret(agg::DataFrame, outdir::String)::String
    "regret_vs_best_static" in names(agg) || return ""
    df = agg[
        (coalesce.(agg.phase_id, "") .== "B6") .&
        (in.(agg.mode, Ref(PPB_ADAPTIVE_MODES))) .&
        .!ismissing.(agg.regret_vs_best_static),
        :,
    ]
    nrow(df) == 0 && return ""

    cases = sort(unique(df.case))
    modes = [m for m in ["outer_inner_adaptive", "full_smart"] if m in unique(df.mode)]
    n_c, n_m = length(cases), length(modes)
    n_c == 0 || n_m == 0 && return ""

    # Worst-case (max) regret per case/mode across whatever thread/process/mc
    # points were tested for it -- one bar per case/mode, not one per point.
    y = [begin
        sub = df[(df.case .== c) .& (df.mode .== m), :]
        nrow(sub) > 0 ? 100.0 * maximum(sub.regret_vs_best_static) : NaN
    end for c in cases, m in modes]

    p = Plots.bar(
        cases, y;
        label         = permutedims(modes),
        ylabel        = "Router regret vs. best static route (%)",
        title         = "B6 — Router Regret Across the Expanded Workload Matrix",
        xrotation     = 25,
        tickfont      = Plots.font(8),
        guidefont     = Plots.font(11),
        titlefont     = Plots.font(11),
        legendfont    = Plots.font(9),
        legend        = :topright,
        size          = (900, 550),
        left_margin   = 14Plots.PlotMeasures.mm,
        bottom_margin = 32Plots.PlotMeasures.mm,
        right_margin  = 24Plots.PlotMeasures.mm,
    )
    Plots.hline!(p, [0.0]; label="", linestyle=:dash, color=:grey, linewidth=1)

    path = joinpath(outdir, "router_regret.png")
    Plots.savefig(p, path)
    return path
end

# ── Top-level plot dispatcher ─────────────────────────────────────────────────

# Runs a single plot-producing closure, isolating its failure so one broken
# plot (bad filter, missing column, empty group, ...) can't prevent the
# others from being generated.
function _ppb_safe_plot(f::Function, label::String)::String
    try
        return f()
    catch err
        @warn "Plot generation failed: $(label)" exception=(err, catch_backtrace())
        return ""
    end
end

function _ppb_write_plots(outdir::String, agg::DataFrame)::Vector{String}
    ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

    paths = String[
        _ppb_safe_plot("constellation_scaling_gravity") do
            _ppb_plot_constellation_scaling(
                agg, "B1",
                "B1 — Constellation Scaling (Gravity Only)",
                "constellation_scaling_gravity.png", outdir)
        end,
        _ppb_safe_plot("constellation_scaling_harmonics") do
            _ppb_plot_constellation_scaling(
                agg, "B2",
                "B2 — Constellation Scaling (L20 Harmonics)",
                "constellation_scaling_harmonics.png", outdir)
        end,
        _ppb_safe_plot("thread_scaling_outer_threads") do
            _ppb_plot_thread_scaling(agg, outdir)
        end,
        _ppb_safe_plot("atmosphere_runtime_comparison") do
            _ppb_plot_atmosphere_runtime(agg, outdir)
        end,
        _ppb_safe_plot("mc_throughput_scaling") do
            _ppb_plot_mc_throughput(agg, outdir)
        end,
        _ppb_safe_plot("profile_comparison") do
            _ppb_plot_profile_comparison(agg, outdir)
        end,
        _ppb_safe_plot("router_regret") do
            _ppb_plot_router_regret(agg, outdir)
        end,
        _ppb_safe_plot("b7_heavy_thread_scaling") do
            _ppb_plot_heavy_thread_scaling(agg, outdir)
        end,
        _ppb_safe_plot("b8_heavy_mc_throughput") do
            _ppb_plot_heavy_mc_throughput(agg, outdir)
        end,
    ]

    filter!(!isempty, paths)
    return paths
end

# ── Report ────────────────────────────────────────────────────────────────────

# Runs a single report section, isolating its failure so a bug in one table
# (bad filter, missing column, ...) can't truncate everything written after
# it — the report file always closes with the sections that did succeed.
function _ppb_safe_section(f::Function, io::IO, label::String)
    try
        f()
    catch err
        @warn "Report section failed: $(label)" exception=(err, catch_backtrace())
        println(io, "_Error generating this section ($(label)); see log for details._")
        println(io)
    end
end

function _ppb_write_report(
    outdir::String,
    agg::DataFrame,
    active_phases::Vector{PPBPhase},
    phase_results::Vector{NamedTuple},
    stamp::String,
)::String
    path = joinpath(outdir, "paper_benchmarks_report_$(stamp).md")
    hw   = ppc_hardware_snapshot()
    open(path, "w") do io
        println(io, "# SpaceAGORA Paper Parallelization Benchmark Report")
        println(io)
        println(io, "Generated: $(now(UTC))")
        println(io, "Machine:   $(hw.machine)")
        println(io, "Julia:     $(hw.julia_version)")
        println(io, "CPU threads available: $(hw.cpu_threads)")
        println(io, "Git commit: $(hw.git_commit)")
        println(io)

        println(io, "## Phase Summary")
        println(io)
        println(io, "| Phase | Label | Status | Runs | Elapsed |")
        println(io, "|-------|-------|--------|------|---------|")
        for r in phase_results
            println(io, "| $(r.id) | $(r.label) | $(r.status) | $(r.runs) | $(_ppb_fmt_elapsed(r.elapsed)) |")
        end
        println(io)

        if nrow(agg) > 0
            _ppb_safe_section(io, "Top Speedups by Case") do
                println(io, "## Top Speedups by Case")
                println(io)
                println(io, "Best speedup achieved across all modes and thread counts for each case,")
                println(io, "excluding the serial baseline.")
                println(io)
                println(io, "| Phase | Case | N sat | Best mode | Threads | Process workers | Speedup |")
                println(io, "|-------|------|------:|-----------|--------:|----------------:|--------:|")

                parallel = agg[agg.mode .!= "serial", :]
                has_pw   = "process_workers" in names(parallel)
                for phase_id in sort(unique(coalesce.(parallel.phase_id, "")))
                    sub = parallel[coalesce.(parallel.phase_id, "") .== phase_id, :]
                    nrow(sub) == 0 && continue
                    for case in sort(unique(sub.case))
                        rows = sub[sub.case .== case, :]
                        nrow(rows) == 0 && continue
                        best_idx = argmax(coalesce.(rows.speedup, -Inf))
                        r = rows[best_idx, :]
                        pw_str   = has_pw ? string(coalesce(r.process_workers, "—")) : "—"
                        sp_str   = ismissing(r.speedup) ? "—" : "$(round(r.speedup; digits=1))×"
                        println(io, "| $(phase_id) | $(r.case) | $(r.n_sat) | $(r.mode) | $(r.thread_count) | $(pw_str) | $(sp_str) |")
                    end
                end
                println(io)
            end

            _ppb_safe_section(io, "Monte Carlo Scaling (B4)") do
                println(io, "## Monte Carlo Scaling (B4)")
                println(io)
                mc = agg[(coalesce.(agg.phase_id, "") .== "B4") .& (agg.mode .== "outer_process"), :]
                if nrow(mc) > 0
                    println(io, "| Case | MC samples | Process workers | Speedup | Process efficiency |")
                    println(io, "|------|:----------:|----------------:|--------:|-------------------:|")
                    mc_sorted = sort(mc, [:case, :mc_samples, :process_workers])
                    for r in eachrow(mc_sorted)
                        pw  = coalesce(r.process_workers, "—")
                        sp  = ismissing(r.speedup) ? "—" : "$(round(r.speedup; digits=1))×"
                        eff = "process_efficiency" in names(mc) && !ismissing(r.process_efficiency) ?
                              "$(round(100 * r.process_efficiency; digits=0))%" : "—"
                        println(io, "| $(r.case) | $(r.mc_samples) | $(pw) | $(sp) | $(eff) |")
                    end
                else
                    println(io, "_No B4 outer_process data._")
                end
                println(io)
            end

            _ppb_safe_section(io, "Router Regret (B6)") do
                println(io, "## Router Regret (B6)")
                println(io)
                println(io, "Regret of the adaptive routing profiles (`outer_inner_adaptive` = R4,")
                println(io, "`full_smart` = R5) relative to the best *static* route (R0-R3) measured on")
                println(io, "the same workload point (case, thread count, process workers, MC samples):")
                println(io)
                println(io, "```")
                println(io, "regret = (T_selected - T_best_static) / T_best_static")
                println(io, "```")
                println(io)
                println(io, "Worst-case (max) regret per case/mode across every point B6 tested for it.")
                println(io, "Negative values mean the adaptive route beat every static route tested; the")
                println(io, "manuscript's own finding that R5 lags ~25% behind the best static route on")
                println(io, "the GRAM workload is exactly what a large positive value here would confirm")
                println(io, "or, if it no longer reproduces, contradict.")
                println(io)
                regret_df = agg[
                    (coalesce.(agg.phase_id, "") .== "B6") .&
                    (in.(agg.mode, Ref(PPB_ADAPTIVE_MODES))) .&
                    .!ismissing.(agg.regret_vs_best_static),
                    :,
                ]
                if nrow(regret_df) > 0
                    println(io, "| Case | Mode | Worst regret | Threads | Process workers | MC samples |")
                    println(io, "|------|------|-------------:|--------:|-----------------:|-----------:|")
                    for case in sort(unique(regret_df.case)), mode in ["outer_inner_adaptive", "full_smart"]
                        sub = regret_df[(regret_df.case .== case) .& (regret_df.mode .== mode), :]
                        nrow(sub) == 0 && continue
                        worst_idx = argmax(sub.regret_vs_best_static)
                        r = sub[worst_idx, :]
                        pw = "process_workers" in names(sub) ? string(coalesce(r.process_workers, "—")) : "—"
                        println(io, "| $(case) | $(mode) | $(round(100 * r.regret_vs_best_static; digits=1))% | $(r.thread_count) | $(pw) | $(r.mc_samples) |")
                    end
                else
                    println(io, "_No B6 router-regret data._")
                end
                println(io)
            end
        end

        println(io, "## Output Files")
        println(io)
        println(io, "Raw and aggregated CSV files are in the same directory as this report.")
        println(io, "Plots: `constellation_scaling_gravity.png`, `constellation_scaling_harmonics.png`,")
        println(io, "`thread_scaling_outer_threads.png`, `atmosphere_runtime_comparison.png`,")
        println(io, "`mc_throughput_scaling.png`, `profile_comparison.png`, `router_regret.png`.")
    end
    return path
end
