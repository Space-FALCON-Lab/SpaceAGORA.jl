const _HYBRID_SCALE_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _HYBRID_SCALE_RUNTIME_SCRIPT = joinpath(@__DIR__, "performance_runtime_analysis.jl")
const _HYBRID_SCALE_DEFAULT_OUTDIR = joinpath(
    _HYBRID_SCALE_REPO_ROOT,
    "output",
    "performance",
    "hybrid_scaling_64_aero_gram",
)
const _HYBRID_SCALE_DEFAULT_CASE = "multi_64_aero_gram"
const _HYBRID_SCALE_DEFAULT_COUNTS = [1, 2, 4, 8, 16, 32, 64]

using CSV
using DataFrames
using Dates
using Statistics

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots

const _HS_PLOT_THEME_APPLIED = Ref(false)

@inline function _hs_env(name::String, default::String)::String
    return strip(get(ENV, name, default))
end

@inline function _hs_bool(name::String, default::Bool)::Bool
    raw = lowercase(_hs_env(name, default ? "1" : "0"))
    return raw in ("1", "true", "yes", "on")
end

function _hs_int_env(name::String, default::Int)::Int
    raw = _hs_env(name, string(default))
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$(name) must be an integer, got '$(raw)'."))
    end
    return value
end

function _hs_counts_env(name::String, default::Vector{Int})::Vector{Int}
    raw = _hs_env(name, join(default, ","))
    isempty(raw) && return copy(default)
    counts = Int[]
    for part in split(raw, ",")
        token = strip(part)
        isempty(token) && continue
        value = try
            parse(Int, token)
        catch
            throw(ArgumentError("$(name) must be a comma-separated list of integers, got '$(raw)'."))
        end
        value > 0 || throw(ArgumentError("$(name) entries must be positive, got $(value)."))
        push!(counts, value)
    end
    isempty(counts) && throw(ArgumentError("$(name) did not contain any positive integers."))
    return sort(unique(counts))
end

function _hs_project()::String
    for candidate in (
        joinpath(_HYBRID_SCALE_REPO_ROOT, ".AGORA"),
        _HYBRID_SCALE_REPO_ROOT,
    )
        if isdir(candidate) && isfile(joinpath(candidate, "Project.toml"))
            return candidate
        end
    end
    return _HYBRID_SCALE_REPO_ROOT
end

function _hs_sysimage()::Union{Nothing, String}
    raw = _hs_env("SPACEAGORA_HYBRID_SYSIMAGE", _hs_env(
        "SPACEAGORA_THREAD_SCALE_SYSIMAGE",
        joinpath(_HYBRID_SCALE_REPO_ROOT, "SpaceAGORA.so"),
    ))
    isempty(raw) && return nothing
    path = abspath(raw)
    return isfile(path) ? path : nothing
end

function _hs_apply_plot_theme!()
    _HS_PLOT_THEME_APPLIED[] && return nothing
    Plots.theme(:ggplot2)
    Plots.default(
        dpi=220,
        lw=2,
        ms=6,
        markerstrokewidth=1.3,
        markerstrokecolor=:black,
        titlefont=Plots.font(20),
        guidefont=Plots.font(14),
        tickfont=Plots.font(11),
        legend_font=Plots.font(10),
        legend_background_color=:white,
        legend_foreground_color=:black,
        gridalpha=0.24,
        framestyle=:box,
    )
    _HS_PLOT_THEME_APPLIED[] = true
    return nothing
end

function _hs_find_latest_raw_csv(outdir::String, profile::String)::Union{Nothing, String}
    isdir(outdir) || return nothing
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

function _hs_command(inner_threads::Int, project::String, profile::String)::Cmd
    julia_bin = Base.julia_cmd().exec[1]
    cmd_parts = String[julia_bin, "--threads=$(inner_threads)"]
    sysimage = _hs_sysimage()
    if sysimage !== nothing
        push!(cmd_parts, "--sysimage=$(sysimage)")
    end
    append!(cmd_parts, String[
        "--project=$(project)",
        _HYBRID_SCALE_RUNTIME_SCRIPT,
        profile,
    ])
    return Cmd(cmd_parts)
end

function _hs_run_worker(
    combo_outdir::String,
    worker_index::Int,
    outer_workers::Int,
    inner_threads::Int,
    case::String,
    profile::String,
    project::String,
)::NamedTuple
    run_outdir = joinpath(combo_outdir, "worker_$(lpad(worker_index, 3, '0'))")
    mkpath(run_outdir)
    cmd = _hs_command(inner_threads, project, profile)
    env_pairs = [
        "SPACEAGORA_PERF_CASES" => case,
        "SPACEAGORA_PERF_OUTDIR" => run_outdir,
        "SPACEAGORA_PERF_INCLUDE_MONTECARLO" => "0",
        "SPACEAGORA_PERF_INCLUDE_MISSION_TIME_SWEEP" => "0",
        "SPACEAGORA_PERF_EXCLUDE_ENTRY_SCENARIOS" => "1",
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => "none",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(inner_threads),
        "SPACEAGORA_PERF_WARMUP_LOGS" => _hs_env("SPACEAGORA_PERF_WARMUP_LOGS", "1"),
        "SPACEAGORA_PERF_WARMUP_OVERRIDE" => _hs_env("SPACEAGORA_PERF_WARMUP_OVERRIDE", "0"),
        "SPACEAGORA_PERF_REPEATS_OVERRIDE" => _hs_env("SPACEAGORA_PERF_REPEATS_OVERRIDE", "1"),
        "SPACEAGORA_PERF_WARMUP_MISSION_SCALE" => _hs_env("SPACEAGORA_PERF_WARMUP_MISSION_SCALE", "0.1"),
    ]
    started_ns = time_ns()
    ok = true
    err_text = ""
    try
        withenv(env_pairs...) do
            run(cmd)
        end
    catch err
        ok = false
        err_text = sprint(showerror, err)
        @warn "[hybrid-scaling] worker failed" outer_workers inner_threads worker_index err_text
    end
    elapsed_s = (time_ns() - started_ns) / 1e9
    raw_csv = _hs_find_latest_raw_csv(run_outdir, profile)
    return (
        worker_index=worker_index,
        ok=ok,
        err_text=err_text,
        elapsed_s=elapsed_s,
        raw_csv=raw_csv === nothing ? "" : raw_csv,
    )
end

function _hs_run_combo(
    outer_workers::Int,
    inner_threads::Int,
    case::String,
    profile::String,
    project::String,
    base_outdir::String,
)::Tuple{DataFrame, DataFrame}
    combo_name = "outer$(outer_workers)_inner$(inner_threads)"
    combo_outdir = joinpath(base_outdir, combo_name)
    mkpath(combo_outdir)
    println()
    println("[hybrid-scaling] outer=$(outer_workers) inner=$(inner_threads) slots=$(outer_workers * inner_threads)")
    combo_started_ns = time_ns()
    worker_results = Vector{NamedTuple}(undef, outer_workers)
    @sync for worker_index in 1:outer_workers
        @async begin
            worker_results[worker_index] = _hs_run_worker(
                combo_outdir,
                worker_index,
                outer_workers,
                inner_threads,
                case,
                profile,
                project,
            )
        end
    end
    combo_elapsed_s = (time_ns() - combo_started_ns) / 1e9

    frames = DataFrame[]
    for result in worker_results
        if !isempty(result.raw_csv) && isfile(result.raw_csv)
            df = CSV.read(result.raw_csv, DataFrame)
            df[!, :outer_workers] .= outer_workers
            df[!, :inner_threads] .= inner_threads
            df[!, :cpu_slots] .= outer_workers * inner_threads
            df[!, :hybrid_worker_index] .= result.worker_index
            df[!, :hybrid_combo_elapsed_s] .= combo_elapsed_s
            df[!, :hybrid_worker_elapsed_s] .= result.elapsed_s
            push!(frames, df)
        end
    end

    worker_df = DataFrame(worker_results)
    worker_df[!, :outer_workers] .= outer_workers
    worker_df[!, :inner_threads] .= inner_threads
    worker_df[!, :cpu_slots] .= outer_workers * inner_threads
    worker_df[!, :combo_elapsed_s] .= combo_elapsed_s

    raw_df = isempty(frames) ? DataFrame() : vcat(frames...; cols=:union)
    println("[hybrid-scaling] outer=$(outer_workers) inner=$(inner_threads) finished in $(round(combo_elapsed_s; digits=2)) s, rows=$(nrow(raw_df))")
    return raw_df, worker_df
end

function _hs_plan(outer_counts::Vector{Int}, inner_counts::Vector{Int}, cpu_budget::Int, allow_oversub::Bool)
    plan = NamedTuple[]
    for outer in outer_counts
        for inner in inner_counts
            slots = outer * inner
            if allow_oversub || slots <= cpu_budget
                push!(plan, (outer_workers=outer, inner_threads=inner, cpu_slots=slots))
            end
        end
    end
    sort!(plan; by = x -> (x.cpu_slots, x.outer_workers, x.inner_threads))
    return plan
end

function _hs_summary(raw_df::DataFrame, worker_df::DataFrame)::DataFrame
    rows = NamedTuple[]
    isempty(worker_df) && return DataFrame(rows)
    for key in sort(unique(zip(worker_df.outer_workers, worker_df.inner_threads)); by=x -> (x[1] * x[2], x[1], x[2]))
        outer, inner = key
        worker_sub = filter(r -> r.outer_workers == outer && r.inner_threads == inner, worker_df)
        raw_sub = isempty(raw_df) ? DataFrame() : filter(r -> r.outer_workers == outer && r.inner_threads == inner, raw_df)
        ok = isempty(raw_sub) ? DataFrame() : filter(r -> r.solve_success === true, raw_sub)
        solve_times = isempty(ok) ? Float64[] : Float64.(ok.solve_time_s)
        combo_elapsed_s = isempty(worker_sub) ? NaN : maximum(Float64.(worker_sub.combo_elapsed_s))
        n_success = length(solve_times)
        mean_solve_s = isempty(solve_times) ? missing : mean(solve_times)
        min_solve_s = isempty(solve_times) ? missing : minimum(solve_times)
        max_solve_s = isempty(solve_times) ? missing : maximum(solve_times)
        throughput = (combo_elapsed_s > 0.0 && n_success > 0) ? n_success / combo_elapsed_s : missing
        push!(rows, (
            outer_workers=outer,
            inner_threads=inner,
            cpu_slots=outer * inner,
            samples_success=n_success,
            samples_total=nrow(raw_sub),
            workers_ok=count(identity, Bool.(worker_sub.ok)),
            workers_total=nrow(worker_sub),
            mean_solve_s=mean_solve_s,
            min_solve_s=min_solve_s,
            max_solve_s=max_solve_s,
            combo_elapsed_s=combo_elapsed_s,
            throughput_solves_per_s=throughput,
        ))
    end
    df = DataFrame(rows)
    df[!, :speedup_vs_1x1] = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))
    df[!, :throughput_speedup_vs_1x1] = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))
    df[!, :efficiency_vs_slots_pct] = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))
    baseline_idx = findfirst(r -> r.outer_workers == 1 && r.inner_threads == 1 && !(r.mean_solve_s isa Missing), eachrow(df))
    if baseline_idx !== nothing
        base_solve = Float64(df.mean_solve_s[baseline_idx])
        base_throughput = df.throughput_solves_per_s[baseline_idx]
        for i in 1:nrow(df)
            if !(df.mean_solve_s[i] isa Missing) && Float64(df.mean_solve_s[i]) > 0.0
                speedup = base_solve / Float64(df.mean_solve_s[i])
                df.speedup_vs_1x1[i] = round(speedup; digits=3)
                df.efficiency_vs_slots_pct[i] = round(100.0 * speedup / Float64(df.cpu_slots[i]); digits=1)
            end
            if !(base_throughput isa Missing) && !(df.throughput_solves_per_s[i] isa Missing) && Float64(base_throughput) > 0.0
                df.throughput_speedup_vs_1x1[i] = round(Float64(df.throughput_solves_per_s[i]) / Float64(base_throughput); digits=3)
            end
        end
    end
    return df
end

function _hs_best_rows(summary_df::DataFrame)
    ok = filter(r -> r.samples_success > 0 && !(r.mean_solve_s isa Missing), summary_df)
    isempty(ok) && return (fastest=nothing, efficient=nothing, economical=nothing)
    fastest = ok[argmin(Float64.(ok.mean_solve_s)), :]
    efficient_candidates = filter(r -> !(r.efficiency_vs_slots_pct isa Missing) && Float64(r.efficiency_vs_slots_pct) >= 50.0, ok)
    efficient = isempty(efficient_candidates) ? nothing : efficient_candidates[argmin(Float64.(efficient_candidates.mean_solve_s)), :]
    best_time = Float64(fastest.mean_solve_s)
    near = filter(r -> Float64(r.mean_solve_s) <= 1.05 * best_time, ok)
    sort!(near, [:cpu_slots, :mean_solve_s])
    economical = near[1, :]
    return (fastest=fastest, efficient=efficient, economical=economical)
end

function _hs_heat_matrix(summary_df::DataFrame, outer_counts::Vector{Int}, inner_counts::Vector{Int}, col::Symbol)
    z = fill(NaN, length(outer_counts), length(inner_counts))
    for row in eachrow(summary_df)
        oi = findfirst(==(row.outer_workers), outer_counts)
        ii = findfirst(==(row.inner_threads), inner_counts)
        (oi === nothing || ii === nothing) && continue
        value = row[col]
        if !(value isa Missing)
            z[oi, ii] = Float64(value)
        end
    end
    return z
end

function _hs_generate_plots(outdir::String, stamp::String, summary_df::DataFrame, outer_counts::Vector{Int}, inner_counts::Vector{Int})::Vector{String}
    paths = String[]
    isempty(summary_df) && return paths
    _hs_apply_plot_theme!()
    xticks = (1:length(inner_counts), string.(inner_counts))
    yticks = (1:length(outer_counts), string.(outer_counts))

    specs = [
        (:mean_solve_s, "Mean Solve Time (s)", "hybrid_mean_solve"),
        (:throughput_speedup_vs_1x1, "Throughput Speedup vs 1x1", "hybrid_throughput_speedup"),
        (:efficiency_vs_slots_pct, "Efficiency vs Slots (%)", "hybrid_efficiency"),
    ]
    for (col, title, stem) in specs
        z = _hs_heat_matrix(summary_df, outer_counts, inner_counts, col)
        p = Plots.heatmap(
            1:length(inner_counts),
            1:length(outer_counts),
            z;
            xlabel="Inner Julia threads per solve",
            ylabel="Outer concurrent solves",
            title=title,
            xticks=xticks,
            yticks=yticks,
            colorbar_title=String(col),
            size=(1200, 820),
            left_margin=16 * Plots.mm,
            bottom_margin=14 * Plots.mm,
        )
        path = joinpath(outdir, "$(stem)_$(stamp).png")
        Plots.savefig(p, path)
        push!(paths, path)
    end

    ok = filter(r -> r.samples_success > 0 && !(r.mean_solve_s isa Missing), summary_df)
    if !isempty(ok)
        p = Plots.scatter(
            ok.cpu_slots,
            Float64.(ok.mean_solve_s);
            group=ok.outer_workers,
            xlabel="CPU slots (outer × inner)",
            ylabel="Mean solve time (s)",
            title="Hybrid Frontier",
            label="outer workers",
            xscale=:log2,
            size=(1200, 760),
            left_margin=16 * Plots.mm,
            bottom_margin=14 * Plots.mm,
        )
        path = joinpath(outdir, "hybrid_frontier_$(stamp).png")
        Plots.savefig(p, path)
        push!(paths, path)
    end
    return paths
end

function _hs_fmt(v)
    v isa Missing && return "n/a"
    v === nothing && return "n/a"
    v isa AbstractFloat && return isfinite(v) ? string(round(v; digits=4)) : "n/a"
    return string(v)
end

function _hs_write_report(path::String, summary_df::DataFrame, case::String, profile::String, cpu_budget::Int, plot_paths::Vector{String})
    best = _hs_best_rows(summary_df)
    open(path, "w") do io
        println(io, "# Hybrid Scaling Study - $(case)")
        println(io)
        println(io, "- Generated (UTC): `$(now(UTC))`")
        println(io, "- Profile: `$(profile)`")
        println(io, "- CPU budget: `$(cpu_budget)`")
        println(io, "- Sysimage: `$(_hs_sysimage() === nothing ? "none" : _hs_sysimage())`")
        println(io, "- Warmup override: `$(get(ENV, "SPACEAGORA_PERF_WARMUP_OVERRIDE", "0"))`")
        println(io, "- Repeats override: `$(get(ENV, "SPACEAGORA_PERF_REPEATS_OVERRIDE", "1"))`")
        println(io, "- Warmup mission scale: `$(get(ENV, "SPACEAGORA_PERF_WARMUP_MISSION_SCALE", "0.1"))`")
        println(io)
        println(io, "## Recommendation")
        println(io)
        if best.fastest === nothing
            println(io, "No successful hybrid samples were collected.")
        else
            f = best.fastest
            e = best.economical
            println(io, "- Fastest: `outer=$(f.outer_workers), inner=$(f.inner_threads), slots=$(f.cpu_slots)`, mean solve `$(_hs_fmt(f.mean_solve_s))` s.")
            println(io, "- Lowest slots within 5% of fastest: `outer=$(e.outer_workers), inner=$(e.inner_threads), slots=$(e.cpu_slots)`, mean solve `$(_hs_fmt(e.mean_solve_s))` s.")
            if best.efficient !== nothing
                eff = best.efficient
                println(io, "- Fastest with at least 50% slot efficiency: `outer=$(eff.outer_workers), inner=$(eff.inner_threads), slots=$(eff.cpu_slots)`, efficiency `$(_hs_fmt(eff.efficiency_vs_slots_pct))`%.")
            end
        end
        println(io)
        println(io, "## Results")
        println(io)
        cols = [
            :outer_workers, :inner_threads, :cpu_slots, :samples_success, :samples_total,
            :workers_ok, :workers_total, :mean_solve_s, :combo_elapsed_s,
            :throughput_solves_per_s, :speedup_vs_1x1, :throughput_speedup_vs_1x1,
            :efficiency_vs_slots_pct,
        ]
        present = [c for c in cols if c in Symbol.(names(summary_df))]
        println(io, "| ", join(String.(present), " | "), " |")
        println(io, "|", join(fill("---:", length(present)), "|"), "|")
        for row in eachrow(summary_df)
            vals = [_hs_fmt(row[col]) for col in present]
            println(io, "| ", join(vals, " | "), " |")
        end
        if !isempty(plot_paths)
            println(io)
            println(io, "## Plots")
            println(io)
            for p in plot_paths
                println(io, "- `$(p)`")
            end
        end
    end
    return nothing
end

function main_hybrid_scaling()
    case = _hs_env("SPACEAGORA_HYBRID_CASE", _HYBRID_SCALE_DEFAULT_CASE)
    profile = _hs_env("SPACEAGORA_HYBRID_PROFILE", "full")
    outdir = _hs_env("SPACEAGORA_HYBRID_OUTDIR", _HYBRID_SCALE_DEFAULT_OUTDIR)
    project = _hs_project()
    cpu_budget = _hs_int_env("SPACEAGORA_HYBRID_CPU_BUDGET", min(Sys.CPU_THREADS, 64))
    cpu_budget > 0 || throw(ArgumentError("SPACEAGORA_HYBRID_CPU_BUDGET must be positive."))
    allow_oversub = _hs_bool("SPACEAGORA_HYBRID_ALLOW_OVERSUBSCRIPTION", false)
    outer_counts = _hs_counts_env("SPACEAGORA_HYBRID_OUTER_COUNTS", _HYBRID_SCALE_DEFAULT_COUNTS)
    inner_counts = _hs_counts_env("SPACEAGORA_HYBRID_INNER_COUNTS", _HYBRID_SCALE_DEFAULT_COUNTS)
    plan = _hs_plan(outer_counts, inner_counts, cpu_budget, allow_oversub)
    isempty(plan) && throw(ArgumentError("No hybrid combinations remain after applying CPU budget $(cpu_budget)."))
    mkpath(outdir)

    println("Hybrid scaling study")
    println("Case:          $(case)")
    println("Profile:       $(profile)")
    println("CPU budget:    $(cpu_budget)")
    println("Combinations:  $(length(plan))")
    println("Sysimage:      $(_hs_sysimage() === nothing ? "none" : _hs_sysimage())")
    println("Output:        $(outdir)")
    println("Project:       $(project)")

    raw_frames = DataFrame[]
    worker_frames = DataFrame[]
    for combo in plan
        raw_df, worker_df = _hs_run_combo(
            combo.outer_workers,
            combo.inner_threads,
            case,
            profile,
            project,
            outdir,
        )
        isempty(raw_df) || push!(raw_frames, raw_df)
        isempty(worker_df) || push!(worker_frames, worker_df)
    end

    all_raw = isempty(raw_frames) ? DataFrame() : vcat(raw_frames...; cols=:union)
    all_workers = isempty(worker_frames) ? DataFrame() : vcat(worker_frames...; cols=:union)
    summary_df = _hs_summary(all_raw, all_workers)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    raw_path = joinpath(outdir, "hybrid_scaling_raw_$(stamp).csv")
    worker_path = joinpath(outdir, "hybrid_scaling_workers_$(stamp).csv")
    summary_path = joinpath(outdir, "hybrid_scaling_summary_$(stamp).csv")
    report_path = joinpath(outdir, "hybrid_scaling_report_$(stamp).md")
    plot_paths = _hs_generate_plots(outdir, stamp, summary_df, outer_counts, inner_counts)

    CSV.write(raw_path, all_raw)
    CSV.write(worker_path, all_workers)
    CSV.write(summary_path, summary_df)
    _hs_write_report(report_path, summary_df, case, profile, cpu_budget, plot_paths)

    best = _hs_best_rows(summary_df)
    println()
    println("Hybrid scaling study complete.")
    println("Raw:     $(raw_path)")
    println("Workers: $(worker_path)")
    println("Summary: $(summary_path)")
    println("Report:  $(report_path)")
    if best.fastest !== nothing
        f = best.fastest
        println("Fastest: outer=$(f.outer_workers), inner=$(f.inner_threads), slots=$(f.cpu_slots), mean_solve=$(_hs_fmt(f.mean_solve_s)) s")
    end
    for p in plot_paths
        println("Plot:    $(p)")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_hybrid_scaling()
end
