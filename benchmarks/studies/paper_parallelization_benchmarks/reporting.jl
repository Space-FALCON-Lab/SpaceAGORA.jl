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
const PPB_ADAPTIVE_MODES = Set(["outer_inner_adaptive", "full_smart", "policy_v2"])

# The phases whose regret figures are the review's answer: B9-B14, the expanded
# per-axis evaluation. B6 is retained as a below-the-floor control (see its
# comment in cli.jl) and is deliberately NOT in this set.
const PPB_ROUTER_PHASES = ["B9", "B10", "B11", "B12", "B13", "B14", "B15"]

# Human-readable name of the workload axis each router phase varies. Used to
# group the per-axis regret summary, which is what point 8 actually asks to see.
const PPB_ROUTER_AXIS_LABELS = Dict(
    "B9"  => "Spacecraft count",
    "B10" => "Atmosphere / GRAM usage",
    "B11" => "Force and actuator model count",
    "B12" => "Interacting vs. independent",
    "B13" => "Thread vs. process budget split",
    "B14" => "Mission duration and output cadence",
    "B15" => "Nested campaign aspect ratio (joint routing)",
    "B6"  => "Small-workload control (below floor)",
)

# Serial wall time below which no routing profile is distinguishable, so a regret
# value computed at that point is measuring fixed per-solve setup rather than the
# router. Derived from the two ends of the existing evidence: B6's whole case list
# sits at 0.017-2.9 s serial and produced 0.0-3.3% regret across every route (i.e.
# nothing), while B7's 8-23 s cases produced -40% to +53%. §7.4 of the findings doc
# puts run-to-run variance on the small measurements at ~16%, which is larger than
# any regret B6 reported. 3 s is the round number between the two regimes; points
# under it are kept and reported but flagged, never quoted as router performance.
const PPB_NOISE_FLOOR_SERIAL_S = 3.0

"""
    _ppb_add_effective_cores!(agg)

Add `effective_cores`: how many cores a measured point could actually occupy.

Regret is only meaningful between points that were given the same machine. The
naive proxy is `thread_count`, and for every thread-backed route that is exactly
right -- the worker subprocess is launched with `--threads=N` and the outer split
runs `Threads.@threads` over it.

It is wrong for the process route. `ppc_run_sample_batch` dispatches that one
across `min(mc_samples, process_workers)` Distributed workers, each its own
`--threads=1` process, and that count does not depend on `thread_count` at all.
So an `outer_process` row recorded at `thread_count=1` really ran on 12 cores.
Scoring a 1-thread adaptive row against it -- which the previous, thread_count-
keyed grouping did -- charged the router for losing a race it entered with a
twelfth of the hardware, and produced a headline "+1428% router regret" on
`independent_1sat_1hr` that was a unit mismatch rather than a measurement.

Serial rows are pinned to 1: `serial` dispatches no outer split and the harness
holds the inner policy off, so its thread_count is a launch parameter it never
spends.
"""
function _ppb_add_effective_cores!(agg::DataFrame)
    nrow(agg) == 0 && return agg
    cols = names(agg)
    ("thread_count" in cols) || return agg
    has_backend = "outer_backend_actual" in cols
    has_procs = "process_workers" in cols
    has_samples = "mc_samples" in cols
    agg[!, :effective_cores] = [
        begin
            backend = has_backend ? lowercase(strip(string(coalesce(agg.outer_backend_actual[i], "")))) : ""
            threads = max(1, Int(agg.thread_count[i]))
            if agg.mode[i] == "serial"
                1
            elseif backend == "process"
                workers = has_procs ? max(1, Int(agg.process_workers[i])) : threads
                samples = has_samples ? max(1, Int(agg.mc_samples[i])) : workers
                min(workers, samples)
            else
                threads
            end
        end
        for i in 1:nrow(agg)
    ]
    return agg
end

function _ppb_add_router_regret!(agg::DataFrame)
    agg[!, :regret_vs_best_static] = Vector{Union{Missing, Float64}}(missing, nrow(agg))
    agg[!, :regret_vs_best_static_oracle] = Vector{Union{Missing, Float64}}(missing, nrow(agg))
    agg[!, :best_static_mode] = Vector{Union{Missing, String}}(missing, nrow(agg))
    agg[!, :below_noise_floor] = fill(false, nrow(agg))
    nrow(agg) == 0 && return agg
    cols = names(agg)
    # Grouped on effective_cores, not thread_count: see _ppb_add_effective_cores!.
    # process_workers is deliberately NOT a key -- it is a pool ceiling that every
    # row in a phase shares, and keying on it would split the process route into
    # its own group with no thread-backed route to compare against.
    budget_col = "effective_cores" in cols ? :effective_cores : :thread_count
    point_key = [k for k in [:phase_id, :case, budget_col, :mc_samples] if string(k) in cols]
    isempty(point_key) && return agg

    # Matched-budget regret: best static route at the SAME thread/worker/sample
    # point. Answers "given this budget, did the router pick the right route?"
    for sub in groupby(agg, point_key)
        static_idx = [
            i for i in 1:nrow(sub)
            if sub.mode[i] in PPB_STATIC_MODES &&
               !ismissing(sub.wall_time_median_s[i]) && sub.wall_time_median_s[i] > 0.0
        ]
        isempty(static_idx) && continue
        best_i = static_idx[argmin(sub.wall_time_median_s[static_idx])]
        best_static = sub.wall_time_median_s[best_i]
        for i in 1:nrow(sub)
            t = sub.wall_time_median_s[i]
            (ismissing(t) || t <= 0.0) && continue
            sub.regret_vs_best_static[i] = (t - best_static) / best_static
            sub.best_static_mode[i] = sub.mode[best_i]
        end
    end

    # Oracle regret: best static route at any budget the router could itself have
    # chosen -- i.e. any thread count *at or below* this row's.
    #
    # The matched-budget figure above holds the width fixed and asks only which
    # route was picked. But R4/R5 also choose how wide to go, so scoring them only
    # against same-width static routes credits them for a width someone else
    # chose: an adaptive route that spreads a workload saturating at 4 threads
    # across 64 shows zero matched-budget regret while wasting 60 cores.
    #
    # The "at or below" bound is what makes this fair rather than merely harsher.
    # `thread_count` is the worker subprocess's Julia --threads value, so a route
    # can always decline to use threads it has but can never exceed them; scoring
    # a 1-thread run against a 12-thread static best would charge it for a budget
    # it was never given (measured on the B7 data, that mistake reported 333%
    # regret for an adaptive route at 1 thread that was in fact within 0.5% of the
    # best route available to it). Under this definition oracle regret collapses
    # to matched-budget regret at the bottom of the ladder, where there is no
    # width choice to make, and diverges at the top, where there is.
    oracle_key = [k for k in [:phase_id, :case, :mc_samples] if string(k) in cols]
    if !isempty(oracle_key) && String(budget_col) in cols
        for sub in groupby(agg, oracle_key)
            for i in 1:nrow(sub)
                t = sub.wall_time_median_s[i]
                (ismissing(t) || t <= 0.0) && continue
                affordable = [
                    sub.wall_time_median_s[j] for j in 1:nrow(sub)
                    if sub.mode[j] in PPB_STATIC_MODES &&
                       !ismissing(sub.wall_time_median_s[j]) && sub.wall_time_median_s[j] > 0.0 &&
                       sub[j, budget_col] <= sub[i, budget_col]
                ]
                isempty(affordable) && continue
                best = minimum(affordable)
                sub.regret_vs_best_static_oracle[i] = (t - best) / best
            end
        end
    end

    # Noise floor, with a fallback for phases that run no serial mode.
    #
    # The floor asks whether a workload is large enough for routes to be
    # distinguishable, and the serial baseline is the natural scale for that. But
    # a phase whose mode ladder omits `serial` -- B13, which compares process and
    # thread splits of a fixed budget and has no serial arm -- gets `missing`
    # there, and treating missing as "below the floor" flagged every one of its
    # points despite wall times of 8.9 to 170 s, the largest in the evaluation.
    # That silently removed a whole axis from the summary.
    #
    # When no serial baseline exists, the slowest configuration measured for that
    # workload stands in for it: if even the worst route resolves faster than the
    # floor, no route can distinguish itself and the point is genuinely too small;
    # if it does not, the workload is measurable whatever the fastest route does.
    # Conservative in the right direction -- it can only admit points, never
    # exclude a real one.
    if "serial_median_s" in names(agg)
        scale_key = [k for k in [:phase_id, :case, :mc_samples] if string(k) in cols]
        fallback = Dict{Any, Float64}()
        if !isempty(scale_key)
            for sub in groupby(agg, scale_key)
                vals = [w for w in sub.wall_time_median_s if !ismissing(w) && w > 0.0]
                isempty(vals) && continue
                fallback[Tuple(sub[1, k] for k in scale_key)] = maximum(vals)
            end
        end
        agg[!, :below_noise_floor] = [
            begin
                s = agg.serial_median_s[i]
                scale = if !ismissing(s)
                    s
                elseif !isempty(scale_key)
                    get(fallback, Tuple(agg[i, k] for k in scale_key), 0.0)
                else
                    0.0
                end
                scale < PPB_NOISE_FLOOR_SERIAL_S
            end
            for i in 1:nrow(agg)
        ]
    end
    return agg
end

# Per-axis regret summary: for each router phase, across every workload point it
# tested that clears the noise floor, the worst and median regret and the fraction
# of points where the adaptive route stayed within 10% of the best static route.
#
# This is the table point 8 is asking for. A worst-case number alone is not a fair
# summary of an adaptive policy (one pathological workload dominates it) and a
# median alone hides exactly the failure the review is pointing at, so all three
# are reported together.
function _ppb_router_regret_summary(agg::DataFrame)::DataFrame
    out = DataFrame(
        phase_id=String[], axis=String[], mode=String[], regret_column=String[],
        n_points=Int[], n_below_floor=Int[], worst_regret=Float64[],
        median_regret=Float64[], frac_within_10pct=Float64[],
    )
    ("regret_vs_best_static" in names(agg) && nrow(agg) > 0) || return out

    for phase in vcat(PPB_ROUTER_PHASES, "B6"), mode in sort(collect(PPB_ADAPTIVE_MODES)),
        (col, colname) in ((:regret_vs_best_static, "matched_budget"),
                           (:regret_vs_best_static_oracle, "oracle"))
        string(col) in names(agg) || continue
        rows = agg[(coalesce.(agg.phase_id, "") .== phase) .& (agg.mode .== mode), :]
        nrow(rows) == 0 && continue
        n_below = count(rows.below_noise_floor)
        scored = rows[.!rows.below_noise_floor .& .!ismissing.(rows[!, col]), :]
        nrow(scored) == 0 && continue
        r = collect(skipmissing(scored[!, col]))
        push!(out, (
            phase, get(PPB_ROUTER_AXIS_LABELS, phase, phase), mode, colname,
            length(r), n_below, maximum(r), median(r), count(<=(0.10), r) / length(r),
        ))
    end
    return out
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

    # The backend the point actually ran on, carried through so regret can be
    # scored against the cores a route really consumed rather than against its
    # nominal thread count. Modes with backend="auto" (R4/R5) resolve this at
    # run time, so it cannot be read off the mode name.
    if "outer_backend_actual" in cols
        backend_df = combine(
            groupby(df, group_keys),
            :outer_backend_actual => (x -> begin
                vals = collect(skipmissing(x))
                isempty(vals) ? "unknown" : String(first(vals))
            end) => :outer_backend_actual,
        )
        agg = leftjoin(agg, backend_df; on=group_keys)
    end
    _ppb_add_effective_cores!(agg)

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
    mode_order = ["serial", "outer_threads", "inner_only", "outer_inner_static", "outer_inner_adaptive", "full_smart", "policy_v2"]
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

# One panel per regret definition (matched-budget and oracle), bars grouped by
# case. Only points that clear the noise floor are plotted: a bar built from
# sub-second workloads is a picture of per-solve setup cost, not of the router,
# and putting it on the same axis as the real measurements invites exactly the
# misreading B6's numbers already caused once.
function _ppb_plot_router_regret(agg::DataFrame, outdir::String)::String
    "regret_vs_best_static" in names(agg) || return ""
    base = agg[
        (in.(coalesce.(agg.phase_id, ""), Ref(Set(PPB_ROUTER_PHASES)))) .&
        (in.(agg.mode, Ref(PPB_ADAPTIVE_MODES))) .&
        .!agg.below_noise_floor,
        :,
    ]
    nrow(base) == 0 && return ""

    modes = [m for m in ["outer_inner_adaptive", "full_smart", "policy_v2"] if m in unique(base.mode)]
    isempty(modes) && return ""

    panels = Any[]
    for (col, title) in ((:regret_vs_best_static, "Matched-budget regret (same thread/worker count)"),
                         (:regret_vs_best_static_oracle, "Oracle regret (best static route at any budget)"))
        string(col) in names(base) || continue
        df = base[.!ismissing.(base[!, col]), :]
        nrow(df) == 0 && continue
        cases = sort(unique(df.case))
        y = [begin
            sub = df[(df.case .== c) .& (df.mode .== m), :]
            nrow(sub) > 0 ? 100.0 * maximum(skipmissing(sub[!, col])) : NaN
        end for c in cases, m in modes]
        pnl = Plots.bar(
            cases, y;
            label      = permutedims(modes),
            ylabel     = "Regret (%)",
            title      = title,
            xrotation  = 30,
            tickfont   = Plots.font(7),
            guidefont  = Plots.font(10),
            titlefont  = Plots.font(10),
            legendfont = Plots.font(8),
            legend     = :topright,
        )
        # 10% is the threshold the per-axis summary reports "fraction within" against.
        Plots.hline!(pnl, [0.0];  label="", linestyle=:dash,  color=:grey, linewidth=1)
        Plots.hline!(pnl, [10.0]; label="", linestyle=:dot,   color=:red,  linewidth=1)
        push!(panels, pnl)
    end
    isempty(panels) && return ""

    p = Plots.plot(
        panels...;
        layout        = (length(panels), 1),
        size          = (1000, 380 * length(panels)),
        left_margin   = 14Plots.PlotMeasures.mm,
        bottom_margin = 34Plots.PlotMeasures.mm,
        right_margin  = 24Plots.PlotMeasures.mm,
        plot_title    = "Router Regret Across the Expanded Workload Matrix (B9-B14)",
        plot_titlefont = Plots.font(12),
    )
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

            _ppb_safe_section(io, "Router Regret by Axis (B9-B14)") do
                println(io, "## Router Regret by Axis (B9-B14)")
                println(io)
                println(io, "Regret of the adaptive routing profiles (`outer_inner_adaptive` = R4,")
                println(io, "`full_smart` = R5) relative to the best *static* route (R0-R3, plus R1b")
                println(io, "where the workload is independent):")
                println(io)
                println(io, "```")
                println(io, "regret = (T_selected - T_best_static) / T_best_static")
                println(io, "```")
                println(io)
                println(io, "Two definitions of `T_best_static`, both reported:")
                println(io)
                println(io, "- **matched_budget** — best static route at the *same* thread/worker count.")
                println(io, "  This is the review's formula literally: given this budget, did the router")
                println(io, "  pick the right route?")
                println(io, "- **oracle** — best static route at *any* budget tested for that case. R4/R5")
                println(io, "  choose the parallel width as well as the route, so the matched-budget")
                println(io, "  figure credits them for a budget someone else picked; the oracle charges")
                println(io, "  an adaptive route for running wide on a workload that saturates narrow.")
                println(io)
                println(io, "Points whose serial baseline is under $(PPB_NOISE_FLOOR_SERIAL_S) s are excluded from these")
                println(io, "statistics (counted under `below floor`): at that size fixed per-solve setup")
                println(io, "dominates, no route is distinguishable, and the resulting regret figure is")
                println(io, "smaller than the run-to-run variance. B6 is listed for exactly that contrast.")
                println(io)
                summary_df = _ppb_router_regret_summary(agg)
                if nrow(summary_df) > 0
                    println(io, "| Phase | Axis | Mode | Regret vs. | Points | Below floor | Worst | Median | Within 10% |")
                    println(io, "|-------|------|------|-----------|-------:|------------:|------:|-------:|-----------:|")
                    for r in eachrow(summary_df)
                        println(io, "| $(r.phase_id) | $(r.axis) | $(r.mode) | $(r.regret_column) | ",
                                "$(r.n_points) | $(r.n_below_floor) | ",
                                "$(round(100 * r.worst_regret; digits=1))% | ",
                                "$(round(100 * r.median_regret; digits=1))% | ",
                                "$(round(Int, 100 * r.frac_within_10pct))% |")
                    end
                else
                    println(io, "_No router-regret data above the noise floor._")
                end
                println(io)

                # Per-case worst points, with the route that actually won and the
                # route the adaptive policy selected. Without the latter two columns
                # a bad regret value says only "something went wrong", not what.
                println(io, "### Worst point per case")
                println(io)
                regret_df = agg[
                    (in.(coalesce.(agg.phase_id, ""), Ref(Set(vcat(PPB_ROUTER_PHASES, "B6"))))) .&
                    (in.(agg.mode, Ref(PPB_ADAPTIVE_MODES))) .&
                    .!ismissing.(agg.regret_vs_best_static),
                    :,
                ]
                if nrow(regret_df) > 0
                    has_plan = "rhs_plan_mode" in names(regret_df)
                    println(io, "| Phase | Case | Mode | Matched | Oracle | Best static was | Router chose | Threads | Workers | MC | Floor |")
                    println(io, "|-------|------|------|--------:|-------:|-----------------|--------------|--------:|--------:|---:|-------|")
                    for phase in vcat(PPB_ROUTER_PHASES, "B6"),
                        case in sort(unique(regret_df[coalesce.(regret_df.phase_id, "") .== phase, :case])),
                        mode in ["outer_inner_adaptive", "full_smart", "policy_v2"]
                        sub = regret_df[
                            (coalesce.(regret_df.phase_id, "") .== phase) .&
                            (regret_df.case .== case) .& (regret_df.mode .== mode), :]
                        nrow(sub) == 0 && continue
                        r = sub[argmax(sub.regret_vs_best_static), :]
                        pw = "process_workers" in names(sub) ? string(coalesce(r.process_workers, "—")) : "—"
                        oracle = ismissing(r.regret_vs_best_static_oracle) ? "—" :
                                 "$(round(100 * r.regret_vs_best_static_oracle; digits=1))%"
                        chose = has_plan ?
                            "$(coalesce(r.rhs_plan_mode, "—"))/$(coalesce(r.rhs_plan_allotment, "—"))" : "—"
                        println(io, "| $(phase) | $(case) | $(mode) | ",
                                "$(round(100 * r.regret_vs_best_static; digits=1))% | $(oracle) | ",
                                "$(coalesce(r.best_static_mode, "—")) | $(chose) | ",
                                "$(r.thread_count) | $(pw) | $(r.mc_samples) | ",
                                "$(r.below_noise_floor ? "below" : "ok") |")
                    end
                else
                    println(io, "_No router-regret data._")
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
