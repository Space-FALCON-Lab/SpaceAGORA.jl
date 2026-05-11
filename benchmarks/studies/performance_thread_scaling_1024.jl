const _THREAD_SCALE_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _THREAD_SCALE_RUNTIME_SCRIPT = joinpath(@__DIR__, "performance_runtime_analysis.jl")
const _THREAD_SCALE_DEFAULT_OUTDIR = joinpath(_THREAD_SCALE_REPO_ROOT, "output", "performance", "thread_scaling_1024")
const _THREAD_SCALE_DEFAULT_CASE = "multi_1024_gravity"
const _THREAD_SCALE_THREAD_COUNTS = [1, 2, 4, 8, 16, 32, 64]

using CSV
using DataFrames
using Dates
using Statistics

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots

const _TS_PLOT_THEME_APPLIED = Ref(false)

function _ts_apply_plot_theme!()
    _TS_PLOT_THEME_APPLIED[] && return nothing
    Plots.theme(:ggplot2)
    Plots.default(
        dpi=220,
        lw=2,
        ms=6,
        markerstrokewidth=1.3,
        markerstrokecolor=:black,
        titlefont=Plots.font(22),
        guidefont=Plots.font(15),
        tickfont=Plots.font(12),
        legend_font=Plots.font(11),
        legend_background_color=:white,
        legend_foreground_color=:black,
        gridalpha=0.24,
        framestyle=:box,
    )
    _TS_PLOT_THEME_APPLIED[] = true
    return nothing
end

@inline function _ts_env(name::String, default::String)::String
    return strip(get(ENV, name, default))
end

@inline function _ts_bool(name::String, default::Bool)::Bool
    raw = lowercase(_ts_env(name, default ? "1" : "0"))
    return raw in ("1", "true", "yes", "on")
end

function _ts_project()::String
    for candidate in (
        joinpath(_THREAD_SCALE_REPO_ROOT, ".AGORA"),
        _THREAD_SCALE_REPO_ROOT,
    )
        if isdir(candidate) && isfile(joinpath(candidate, "Project.toml"))
            return candidate
        end
    end
    return _THREAD_SCALE_REPO_ROOT
end

function _ts_find_latest_raw_csv(outdir::String, profile::String)::Union{Nothing, String}
    prefix = "runtime_raw_$(profile)_"
    candidates = [
        joinpath(outdir, f)
        for f in readdir(outdir)
        if startswith(f, prefix) && endswith(f, ".csv")
    ]
    isempty(candidates) && return nothing
    sort!(candidates; by = p -> stat(p).mtime)
    return last(candidates)
end

function _ts_run_one(
    n_threads::Int,
    case::String,
    profile::String,
    project::String,
    base_outdir::String
)::Union{Nothing, DataFrame}
    run_outdir = joinpath(base_outdir, "$(n_threads)t")
    mkpath(run_outdir)

    println()
    println("[thread-scaling] threads=$(n_threads)  case=$(case)  profile=$(profile)")

    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([
        julia_bin,
        "--threads=$(n_threads)",
        "--project=$(project)",
        _THREAD_SCALE_RUNTIME_SCRIPT,
        profile,
    ])

    env_pairs = [
        "SPACEAGORA_PERF_CASES"                   => case,
        "SPACEAGORA_PERF_OUTDIR"                  => run_outdir,
        "SPACEAGORA_PERF_INCLUDE_MONTECARLO"      => "0",
        "SPACEAGORA_PERF_INCLUDE_MISSION_TIME_SWEEP" => "0",
        "SPACEAGORA_PERF_EXCLUDE_ENTRY_SCENARIOS" => "1",
        "SPACEAGORA_PERF_PARALLEL_BACKEND"        => "none",
        "SPACEAGORA_INNER_THREAD_BUDGET"          => "0",
    ]

    started_ns = time_ns()
    try
        withenv(env_pairs...) do
            run(cmd)
        end
    catch err
        @warn "[thread-scaling] subprocess for threads=$(n_threads) failed: $(err)"
        return nothing
    end
    elapsed_s = (time_ns() - started_ns) / 1e9
    println("[thread-scaling] threads=$(n_threads) finished in $(round(elapsed_s; digits=2)) s")

    raw_csv = _ts_find_latest_raw_csv(run_outdir, profile)
    if raw_csv === nothing
        @warn "[thread-scaling] no raw CSV found in $(run_outdir) for profile=$(profile)"
        return nothing
    end
    return CSV.read(raw_csv, DataFrame)
end

function _ts_build_summary(all_raw::DataFrame)::DataFrame
    rows = NamedTuple[]
    thread_col = :julia_threads in Symbol.(names(all_raw)) ? :julia_threads : nothing

    thread_col === nothing && @warn "[thread-scaling] 'julia_threads' column not found in raw data; summary may be incomplete."

    groups = thread_col === nothing ? [(missing, all_raw)] :
        [(t, filter(r -> r[thread_col] == t, all_raw)) for t in sort(unique(all_raw[!, thread_col]))]

    for (n_threads, group) in groups
        for scenario in sort(unique(group.scenario))
            sub = filter(r -> r.scenario == scenario, group)
            ok = filter(r -> r.solve_success, sub)
            n_ok = nrow(ok)
            n_total = nrow(sub)

            if n_ok == 0
                push!(rows, (
                    n_threads          = n_threads,
                    scenario           = scenario,
                    samples_success    = 0,
                    samples_total      = n_total,
                    mean_solve_s       = missing,
                    min_solve_s        = missing,
                    max_solve_s        = missing,
                    solver_sequence    = missing,
                    fallback_used      = missing,
                    mean_solve_allocs  = missing,
                    speedup_vs_1t      = missing,
                    parallel_efficiency_pct = missing,
                ))
                continue
            end

            times = Float64.(ok.solve_time_s)
            solver_seq = :solver_sequence in names(ok) ? first(ok.solver_sequence) : missing
            fallback = :solver_fallback_used in names(ok) ? any(skipmissing(Bool.(ok.solver_fallback_used))) : missing
            allocs = :solve_alloc_calls in names(ok) ? mean(Float64.(ok.solve_alloc_calls)) : missing

            push!(rows, (
                n_threads          = n_threads,
                scenario           = scenario,
                samples_success    = n_ok,
                samples_total      = n_total,
                mean_solve_s       = mean(times),
                min_solve_s        = minimum(times),
                max_solve_s        = maximum(times),
                solver_sequence    = solver_seq,
                fallback_used      = fallback,
                mean_solve_allocs  = allocs,
                speedup_vs_1t      = missing,
                parallel_efficiency_pct = missing,
            ))
        end
    end

    df = DataFrame(rows)
    df[!, :speedup_vs_1t] = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))
    df[!, :parallel_efficiency_pct] = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))

    # Fill speedup and efficiency relative to the 1-thread row.
    for scenario in unique(df.scenario)
        mask_1t = (df.scenario .== scenario) .& (df.n_threads .=== 1)
        any(mask_1t) || continue
        baseline = df.mean_solve_s[findfirst(mask_1t)]
        (baseline isa Missing || !isfinite(Float64(baseline)) || Float64(baseline) <= 0.0) && continue
        for i in 1:nrow(df)
            df.scenario[i] == scenario || continue
            nt = df.n_threads[i]
            (nt isa Missing || nt <= 0) && continue
            t = df.mean_solve_s[i]
            (t isa Missing || !isfinite(Float64(t)) || Float64(t) <= 0.0) && continue
            speedup = Float64(baseline) / Float64(t)
            df.speedup_vs_1t[i] = round(speedup; digits=3)
            df.parallel_efficiency_pct[i] = round(100.0 * speedup / Float64(nt); digits=1)
        end
    end

    sort!(df, [:scenario, :n_threads])
    return df
end

function _ts_generate_plots(
    outdir::String,
    stamp::String,
    summary_df::DataFrame,
    thread_counts::Vector{Int},
    case::String,
    profile::String,
)::Vector{String}
    paths = String[]
    _ts_apply_plot_theme!()

    case_df = filter(r -> r.scenario == case, summary_df)
    isempty(case_df) && return paths

    xs = Int[]
    ys_time = Float64[]
    ys_eff = Float64[]

    for n in thread_counts
        row_idx = findfirst(r -> r.n_threads == n, eachrow(case_df))
        row_idx === nothing && continue
        row = case_df[row_idx, :]
        t = row.mean_solve_s
        e = row.parallel_efficiency_pct
        (t isa Missing || !isfinite(Float64(t))) && continue
        push!(xs, n)
        push!(ys_time, Float64(t))
        push!(ys_eff, (e isa Missing || !isfinite(Float64(e))) ? NaN : Float64(e))
    end

    isempty(xs) && return paths

    x_ticks = (xs, string.(xs))
    baseline_time = ys_time[1]
    ideal_time = [baseline_time / n for n in xs]

    # Panel 1: solve time vs threads (log-log).
    p1 = Plots.plot(
        xs, ys_time;
        label="Measured",
        color="#2f7fc1",
        marker=:circle,
        xscale=:log2,
        yscale=:log10,
        xticks=x_ticks,
        title="Solve Time vs Thread Count\n($(case), $(profile) profile)",
        xlabel="Julia threads",
        ylabel="Mean solve time (s)",
        size=(1700, 820),
        left_margin=20 * Plots.mm,
        right_margin=20 * Plots.mm,
        top_margin=10 * Plots.mm,
        bottom_margin=16 * Plots.mm,
        legend=:topright,
        framestyle=:box,
        gridalpha=0.24,
        legend_background_color=:white,
    )
    Plots.plot!(p1, xs, ideal_time;
        label="Ideal (linear scaling)",
        color=:black,
        linestyle=:dash,
        marker=:none,
    )

    # Panel 2: parallel efficiency vs threads (linear).
    eff_finite = replace(ys_eff, NaN => missing)
    p2 = Plots.plot(
        xs, ys_eff;
        label="Efficiency",
        color="#d67c1c",
        marker=:square,
        xticks=x_ticks,
        ylims=(0, max(115.0, maximum(filter(isfinite, ys_eff); init=100.0) + 10.0)),
        title="Parallel Efficiency",
        xlabel="Julia threads",
        ylabel="Efficiency (%)",
        size=(1700, 820),
        left_margin=20 * Plots.mm,
        right_margin=20 * Plots.mm,
        top_margin=10 * Plots.mm,
        bottom_margin=16 * Plots.mm,
        legend=:topright,
        framestyle=:box,
        gridalpha=0.24,
        legend_background_color=:white,
    )
    Plots.hline!(p2, [100.0]; color=:black, linestyle=:dash, label="100% (ideal)")

    combined = Plots.plot(p1, p2;
        layout=(2, 1),
        size=(1700, 1500),
        left_margin=22 * Plots.mm,
        right_margin=22 * Plots.mm,
        top_margin=10 * Plots.mm,
        bottom_margin=18 * Plots.mm,
    )

    combined_path = joinpath(outdir, "thread_scaling_plot_$(profile)_$(stamp).png")
    try
        Plots.savefig(combined, combined_path)
        push!(paths, combined_path)
        println("[thread-scaling] plot saved: $(combined_path)")
    catch err
        @warn "[thread-scaling] failed to save plot: $(err)"
    end

    return paths
end

function _ts_write_report(
    path::String,
    summary_df::DataFrame,
    thread_counts::Vector{Int},
    case::String,
    profile::String;
    plot_paths::Vector{String}=String[],
)
    open(path, "w") do io
        println(io, "# Thread Scaling Study — $(case)")
        println(io)
        println(io, "- Generated (UTC): `$(now(UTC))`")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Profile: `$(profile)`")
        println(io, "- Case: `$(case)`")
        println(io, "- Thread counts measured: $(join(thread_counts, ", "))")
        println(io, "- Host CPU threads available: `$(Sys.CPU_THREADS)`")
        println(io)
        println(io, "## Scaling Table")
        println(io)
        println(io, "- `speedup_vs_1t`: 1-thread mean solve time / N-thread mean solve time.")
        println(io, "- `parallel_efficiency_pct`: `100 × speedup / n_threads`. 100% = perfect linear scaling.")
        println(io)

        col_order = [
            :n_threads, :scenario, :samples_success, :samples_total,
            :mean_solve_s, :min_solve_s, :max_solve_s,
            :speedup_vs_1t, :parallel_efficiency_pct,
            :solver_sequence, :fallback_used, :mean_solve_allocs,
        ]
        present = [c for c in col_order if c in Symbol.(names(summary_df))]
        header = join(String.(present), " | ")
        println(io, "| ", header, " |")
        println(io, "|", join(fill("---:", length(present)), "|"), "|")
        for row in eachrow(summary_df)
            vals = map(present) do col
                v = row[col]
                v isa Missing && return "n/a"
                v isa Bool && return v ? "true" : "false"
                v isa Float64 && return isfinite(v) ? string(round(v; digits=4)) : "n/a"
                string(v)
            end
            println(io, "| ", join(vals, " | "), " |")
        end
        println(io)
        println(io, "## Reproducibility")
        println(io)
        println(io, "Each row is a fresh subprocess:")
        println(io)
        println(io, "```bash")
        for n in thread_counts
            println(io,
                "SPACEAGORA_PERF_CASES=$(case) " *
                "SPACEAGORA_PERF_INCLUDE_MONTECARLO=0 " *
                "SPACEAGORA_PERF_INCLUDE_MISSION_TIME_SWEEP=0 " *
                "SPACEAGORA_PERF_EXCLUDE_ENTRY_SCENARIOS=1 " *
                "SPACEAGORA_PERF_PARALLEL_BACKEND=none " *
                "julia --threads=$(n) --project=. benchmarks/studies/performance_runtime_analysis.jl $(profile)"
            )
        end
        println(io, "```")
        println(io)

        if !isempty(plot_paths)
            println(io, "## Plots")
            println(io)
            for p in plot_paths
                println(io, "- `$(p)`")
            end
            println(io)
        end
    end
end

function main_thread_scaling()
    outdir  = _ts_env("SPACEAGORA_THREAD_SCALE_OUTDIR", _THREAD_SCALE_DEFAULT_OUTDIR)
    profile = _ts_env("SPACEAGORA_THREAD_SCALE_PROFILE", "full")
    case    = _ts_env("SPACEAGORA_THREAD_SCALE_CASE", _THREAD_SCALE_DEFAULT_CASE)
    project = _ts_project()

    thread_counts = _THREAD_SCALE_THREAD_COUNTS
    mkpath(outdir)

    println("Thread scaling study")
    println("Case:          $(case)")
    println("Profile:       $(profile)")
    println("Thread counts: $(join(thread_counts, ", "))")
    println("Output:        $(outdir)")
    println("Project:       $(project)")

    frames = DataFrame[]
    for n in thread_counts
        df = _ts_run_one(n, case, profile, project, outdir)
        df === nothing || push!(frames, df)
    end

    if isempty(frames)
        @warn "[thread-scaling] no results collected."
        return
    end

    all_raw = vcat(frames...; cols=:union)
    summary_df = _ts_build_summary(all_raw)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    raw_path     = joinpath(outdir, "thread_scaling_raw_$(stamp).csv")
    summary_path = joinpath(outdir, "thread_scaling_summary_$(stamp).csv")
    report_path  = joinpath(outdir, "thread_scaling_report_$(stamp).md")

    plot_paths = _ts_generate_plots(outdir, stamp, summary_df, thread_counts, case, profile)

    CSV.write(raw_path, all_raw)
    CSV.write(summary_path, summary_df)
    _ts_write_report(report_path, summary_df, thread_counts, case, profile; plot_paths=plot_paths)

    println()
    println("Thread scaling study complete.")
    println("Combined raw: $(raw_path)")
    println("Summary:      $(summary_path)")
    println("Report:       $(report_path)")
    for p in plot_paths
        println("Plot:         $(p)")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_thread_scaling()
end
