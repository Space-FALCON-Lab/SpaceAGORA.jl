# Monte Carlo outer-parallelism distributed-worker scaling benchmark.
#
# For each worker count N, spawns a coordinator subprocess (performance_mc_thread_scaling_worker.jl)
# which starts N distributed worker processes via addprocs — the same pattern used by
# aerobraking_perturbation_mc/main.jl — and dispatches N independent single-satellite
# simulations to those workers via remotecall_fetch. Wall-clock time is recorded.
# The key metric is:
#
#   time_per_case = wall_time / N
#
# The optimal worker count is the N that minimises time_per_case (maximises throughput).
# This is distinct from the inner-thread scaling benchmarks, which measure how a single
# simulation speeds up with more threads; here all parallelism is outer (independent cases
# across separate processes).
#
# Usage:
#   julia --project=. benchmarks/studies/performance_mc_thread_scaling.jl
#
# Environment variables:
#   SPACEAGORA_MC_THREAD_COUNTS    comma-separated worker counts (default: powers of two
#                                  up to the physical core count, plus the physical and
#                                  logical core counts)
#   SPACEAGORA_MC_THREAD_OUTDIR    output directory
#   SPACEAGORA_MC_THREAD_NORBITS   orbital periods per case (default: 10)
#   SPACEAGORA_MC_THREAD_BATCHES   cases per worker (default: 1; use >=2 for steadier
#                                  throughput estimates that amortise dispatch jitter)
#   SPACEAGORA_MC_PHYSICAL_CORES   override detected physical core count
#   SPACEAGORA_MC_THREAD_SYSIMAGE  path to sysimage .so (auto-detected if omitted)
#   SPACEAGORA_MC_THREAD_SMOKE     set to "1" for a minimal 1-worker smoke run

const _MC_TS_REPO_ROOT    = normpath(joinpath(@__DIR__, "..", ".."))
const _MC_TS_WORKER_SCRIPT = joinpath(@__DIR__, "performance_mc_thread_scaling_worker.jl")
const _MC_TS_DEFAULT_OUTDIR = joinpath(_MC_TS_REPO_ROOT, "output", "performance", "mc_thread_scaling")

# Physical core count, distinct from Sys.CPU_THREADS (logical, includes SMT).
# Returns (count, exact): exact=false means we fell back to the logical count.
function _mc_ts_physical_cores()::Tuple{Int, Bool}
    override = strip(get(ENV, "SPACEAGORA_MC_PHYSICAL_CORES", ""))
    isempty(override) || return (parse(Int, override), true)
    if Sys.isapple()
        try
            return (parse(Int, strip(read(`sysctl -n hw.physicalcpu`, String))), true)
        catch
        end
    elseif Sys.islinux()
        try
            pairs = Set{String}()
            for line in eachline(`lscpu -p=core,socket`)
                startswith(line, "#") && continue
                isempty(strip(line)) && continue
                push!(pairs, strip(line))
            end
            isempty(pairs) || return (length(pairs), true)
        catch
        end
    end
    return (Sys.CPU_THREADS, false)
end

# Default sweep: powers of two up to the physical core count, the physical count
# itself, and the logical count when SMT is present (so the SMT plateau is visible
# without wasting runs far beyond it).
function _mc_ts_default_thread_counts()::Vector{Int}
    physical, _ = _mc_ts_physical_cores()
    counts = Int[]
    n = 1
    while n < physical
        push!(counts, n)
        n *= 2
    end
    push!(counts, physical)
    Sys.CPU_THREADS > physical && push!(counts, Sys.CPU_THREADS)
    return sort!(unique!(counts))
end

using CSV
using DataFrames
using Dates
using Printf
using Statistics
using TOML

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots

const _MC_TS_PLOT_THEME_APPLIED = Ref(false)

function _mc_ts_apply_plot_theme!()
    _MC_TS_PLOT_THEME_APPLIED[] && return nothing
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
    _MC_TS_PLOT_THEME_APPLIED[] = true
    return nothing
end

@inline function _mc_ts_env(name::String, default::String)::String
    return strip(get(ENV, name, default))
end

@inline function _mc_ts_bool(name::String, default::Bool)::Bool
    raw = lowercase(_mc_ts_env(name, default ? "1" : "0"))
    return raw in ("1", "true", "yes", "on")
end

function _mc_ts_thread_counts()::Vector{Int}
    raw = _mc_ts_env("SPACEAGORA_MC_THREAD_COUNTS", join(_mc_ts_default_thread_counts(), ","))
    counts = Int[]
    for part in split(raw, ",")
        t = strip(part)
        isempty(t) && continue
        n = parse(Int, t)
        n > 0 || throw(ArgumentError("thread count must be positive, got $(n)"))
        push!(counts, n)
    end
    isempty(counts) && throw(ArgumentError("SPACEAGORA_MC_THREAD_COUNTS contained no valid entries"))
    return sort!(unique!(counts))
end

function _mc_ts_project()::String
    for candidate in (
        joinpath(_MC_TS_REPO_ROOT, ".AGORA"),
        _MC_TS_REPO_ROOT,
    )
        if isdir(candidate) && isfile(joinpath(candidate, "Project.toml"))
            return candidate
        end
    end
    return _MC_TS_REPO_ROOT
end

function _mc_ts_sysimage()::Union{Nothing, String}
    raw = _mc_ts_env("SPACEAGORA_MC_THREAD_SYSIMAGE",
        _mc_ts_env("SPACEAGORA_THREAD_SCALE_SYSIMAGE",
            joinpath(_MC_TS_REPO_ROOT, "SpaceAGORA.so")))
    isempty(raw) && return nothing
    path = abspath(raw)
    return isfile(path) ? path : nothing
end

function _mc_ts_run_one(
    n_threads::Int,
    n_orbits::Int,
    project::String,
    base_outdir::String,
)::Union{Nothing, NamedTuple}
    run_outdir = joinpath(base_outdir, "$(n_threads)t")
    mkpath(run_outdir)

    println()
    @printf("[mc-thread-scaling] threads=%-3d  norbits=%d\n", n_threads, n_orbits)
    flush(stdout)

    julia_bin = Base.julia_cmd().exec[1]

    # The coordinator process only dispatches work via @async/remotecall_fetch and needs
    # just one thread. Workers are separate processes started by the coordinator via addprocs;
    # their thread count and inner-parallelism env are set there.
    cmd_parts = String[
        julia_bin,
        "--threads=1",
        "--gcthreads=1,1",
    ]
    sysimage = _mc_ts_sysimage()
    if sysimage !== nothing
        push!(cmd_parts, "--sysimage=$(sysimage)")
    end
    append!(cmd_parts, [
        "--project=$(project)",
        _MC_TS_WORKER_SCRIPT,
    ])
    cmd = Cmd(cmd_parts)

    elapsed_s = try
        withenv(
            "SPACEAGORA_MC_THREAD_OUTDIR"  => run_outdir,
            "SPACEAGORA_MC_THREAD_NORBITS" => string(n_orbits),
            "SPACEAGORA_MC_DIST_NWORKERS"  => string(n_threads),
            "SPACEAGORA_MC_THREAD_BATCHES" => _mc_ts_env("SPACEAGORA_MC_THREAD_BATCHES", "1"),
        ) do
            @elapsed run(cmd)
        end
    catch err
        @warn "[mc-thread-scaling] subprocess for threads=$(n_threads) failed: $(err)"
        return nothing
    end

    @printf("[mc-thread-scaling] threads=%-3d  subprocess wall=%.2f s\n", n_threads, elapsed_s)
    flush(stdout)

    result_csv = joinpath(run_outdir, "mc_ts_worker_$(n_threads)t.csv")
    if !isfile(result_csv)
        @warn "[mc-thread-scaling] expected output not found: $(result_csv)"
        return nothing
    end

    df = CSV.read(result_csv, DataFrame)
    nrow(df) == 0 && return nothing
    row = df[1, :]
    row_names = Symbol.(names(df))
    return (
        julia_threads      = Int(row.julia_threads),
        n_concurrent       = Int(row.n_concurrent),
        n_cases            = :n_cases in row_names ? Int(row.n_cases) : Int(row.n_concurrent),
        batches_per_worker = :batches_per_worker in row_names ? Int(row.batches_per_worker) : 1,
        wall_time_s        = Float64(row.wall_time_s),
        time_per_case_s    = Float64(row.time_per_case_s),
        throughput_cases_s = Float64(row.throughput_cases_s),
        mean_case_s        = Float64(row.mean_case_s),
        min_case_s         = Float64(row.min_case_s),
        max_case_s         = Float64(row.max_case_s),
        imbalance_pct      = Float64(row.imbalance_pct),
        mean_gc_pct        = :mean_gc_pct in row_names ? Float64(row.mean_gc_pct) : NaN,
        n_orbits           = Int(row.n_orbits),
        case_id            = String(row.case_id),
    )
end

function _mc_ts_build_summary(rows::Vector)::DataFrame
    df = DataFrame(rows)

    df[!, :speedup_throughput]     = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))
    df[!, :efficiency_pct]         = Vector{Union{Missing, Float64}}(fill(missing, nrow(df)))
    df[!, :is_optimal]             = fill(false, nrow(df))

    # Baseline is 1-thread row (throughput = 1 / time_per_case).
    idx_1t = findfirst(r -> r.julia_threads == 1, eachrow(df))
    if idx_1t !== nothing
        base_tpc = df.time_per_case_s[idx_1t]   # single-threaded time-per-case
        if isfinite(base_tpc) && base_tpc > 0.0
            for i in 1:nrow(df)
                n  = df.julia_threads[i]
                tp = df.time_per_case_s[i]
                (isfinite(tp) && tp > 0.0) || continue
                speedup = base_tpc / tp
                df.speedup_throughput[i] = round(speedup; digits=3)
                df.efficiency_pct[i]     = round(100.0 * speedup / n; digits=1)
            end
        end
    end

    # Optimal = lowest time_per_case.
    valid_tpc = filter(isfinite, df.time_per_case_s)
    if !isempty(valid_tpc)
        best_tpc = minimum(valid_tpc)
        for i in 1:nrow(df)
            isfinite(df.time_per_case_s[i]) || continue
            df.is_optimal[i] = (df.time_per_case_s[i] ≈ best_tpc)
        end
    end

    sort!(df, :julia_threads)
    return df
end

function _mc_ts_optimal_row(summary_df::DataFrame)::Union{Nothing, DataFrameRow}
    opt_rows = filter(r -> r.is_optimal, eachrow(summary_df))
    isempty(opt_rows) && return nothing
    # Among optimal rows (may be a tie due to ≈), prefer the one with highest throughput.
    return first(sort(collect(opt_rows); by=r -> -r.throughput_cases_s))
end

function _mc_ts_generate_plots(
    outdir::String,
    stamp::String,
    summary_df::DataFrame,
)::Vector{String}
    paths = String[]
    isempty(summary_df) && return paths
    _mc_ts_apply_plot_theme!()

    xs     = Int.(summary_df.julia_threads)
    tpc    = Float64.(summary_df.time_per_case_s)
    thput  = Float64.(summary_df.throughput_cases_s)
    eff    = [v isa Missing ? NaN : Float64(v) for v in summary_df.efficiency_pct]

    x_ticks = (xs, string.(xs))

    # Ideal lines anchored at 1-thread result.
    idx_1t = findfirst(==(1), xs)
    ideal_tpc   = idx_1t !== nothing ? [tpc[idx_1t] / n  for n in xs] : fill(NaN, length(xs))
    ideal_thput = idx_1t !== nothing ? [thput[idx_1t] * n for n in xs] : fill(NaN, length(xs))

    # Mark optimal point.
    opt_row = _mc_ts_optimal_row(summary_df)
    opt_n   = opt_row !== nothing ? opt_row.julia_threads : -1

    physical_cores, physical_exact = _mc_ts_physical_cores()
    show_core_line = physical_exact && minimum(xs) <= physical_cores <= maximum(xs)
    core_label = "Physical cores ($(physical_cores))"

    # Panel 1: time-per-case vs workers.
    p1 = Plots.plot(
        xs, tpc;
        label="Measured",
        color="#2f7fc1",
        marker=:circle,
        xscale=:log2,
        yscale=:log10,
        xticks=x_ticks,
        title="MC Time-per-Case vs Distributed Worker Count",
        xlabel="Distributed workers",
        ylabel="Time per case (s)",
        size=(1700, 820),
        left_margin=20Plots.mm,
        right_margin=20Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=16Plots.mm,
        legend=:topright,
    )
    Plots.plot!(p1, xs, ideal_tpc;
        label="Ideal (linear speedup)",
        color=:black,
        linestyle=:dash,
        marker=:none,
    )
    if opt_n > 0
        idx_opt = findfirst(==(opt_n), xs)
        if idx_opt !== nothing
            Plots.scatter!(p1, [xs[idx_opt]], [tpc[idx_opt]];
                label="Optimal (n=$(opt_n))",
                color="#e84545",
                marker=:star5,
                ms=10,
            )
        end
    end
    show_core_line && Plots.vline!(p1, [physical_cores]; color=:gray, linestyle=:dot, label=core_label)

    # Panel 2: throughput (cases/s) vs workers.
    p2 = Plots.plot(
        xs, thput;
        label="Measured",
        color="#2aab5a",
        marker=:diamond,
        xscale=:log2,
        xticks=x_ticks,
        title="MC Throughput vs Distributed Worker Count",
        xlabel="Distributed workers",
        ylabel="Throughput (cases / s)",
        size=(1700, 820),
        left_margin=20Plots.mm,
        right_margin=20Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=16Plots.mm,
        legend=:topleft,
    )
    show_core_line && Plots.vline!(p2, [physical_cores]; color=:gray, linestyle=:dot, label=core_label)
    Plots.plot!(p2, xs, ideal_thput;
        label="Ideal (linear throughput)",
        color=:black,
        linestyle=:dash,
        marker=:none,
    )

    # Panel 3: parallel efficiency.
    eff_finite = filter(isfinite, eff)
    p3 = Plots.plot(
        xs, eff;
        label="Efficiency",
        color="#d67c1c",
        marker=:square,
        xticks=x_ticks,
        ylims=(0, max(115.0, isempty(eff_finite) ? 110.0 : maximum(eff_finite) + 10.0)),
        title="MC Parallel Efficiency (throughput basis)",
        xlabel="Distributed workers",
        ylabel="Efficiency (%)",
        size=(1700, 820),
        left_margin=20Plots.mm,
        right_margin=20Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=16Plots.mm,
        legend=:topright,
    )
    Plots.hline!(p3, [100.0]; color=:black, linestyle=:dash, label="100% (ideal)")
    show_core_line && Plots.vline!(p3, [physical_cores]; color=:gray, linestyle=:dot, label=core_label)

    combined = Plots.plot(p1, p2, p3;
        layout=(3, 1),
        size=(1700, 2200),
        left_margin=22Plots.mm,
        right_margin=22Plots.mm,
        top_margin=10Plots.mm,
        bottom_margin=18Plots.mm,
    )

    plot_path = joinpath(outdir, "mc_thread_scaling_plot_$(stamp).png")
    try
        Plots.savefig(combined, plot_path)
        push!(paths, plot_path)
        println("[mc-thread-scaling] plot saved: $(plot_path)")
    catch err
        @warn "[mc-thread-scaling] failed to save plot: $(err)"
    end

    return paths
end

function _mc_ts_write_report(
    path::String,
    summary_df::DataFrame,
    thread_counts::Vector{Int},
    n_orbits::Int;
    plot_paths::Vector{String}=String[],
)
    opt_row = _mc_ts_optimal_row(summary_df)
    opt_n   = opt_row !== nothing ? opt_row.julia_threads : -1

    open(path, "w") do io
        println(io, "# Monte Carlo Distributed Worker Scaling Study")
        println(io)
        physical_cores, physical_exact = _mc_ts_physical_cores()
        println(io, "- Generated (UTC): `$(now(UTC))`")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Worker counts measured: $(join(thread_counts, ", "))")
        println(io, "- Host physical cores: `$(physical_cores)`$(physical_exact ? "" : " (estimate: fell back to logical count)")")
        println(io, "- Host CPU threads available (logical): `$(Sys.CPU_THREADS)`")
        println(io, "- Workload: $(nrow(summary_df) > 0 ? summary_df.case_id[1] : "unknown")")
        println(io, "- Orbits per case: $(n_orbits)")
        sysimage = _mc_ts_sysimage()
        println(io, "- Sysimage: `$(sysimage === nothing ? "none" : sysimage)`")
        println(io)
        println(io, "## Optimal Thread Count")
        println(io)
        if opt_n > 0 && opt_row !== nothing
            println(io, "**$(opt_n) distributed workers** minimises `time_per_case` at **$(round(opt_row.time_per_case_s; digits=3)) s/case** " *
                "(throughput: $(round(opt_row.throughput_cases_s; digits=3)) cases/s).")
            println(io)
            println(io, "Use this when running MC studies with distributed outer parallelism:")
            println(io)
            println(io, "```bash")
            println(io, "julia --project=. \\")
            println(io, "  benchmarks/studies/aerobraking_perturbation_mc/main.jl --procs $(opt_n)")
            println(io, "```")
            println(io)
            println(io, "The optimal count has been written to `mc_optimal_threads.toml` in the output directory.")
        else
            println(io, "Could not determine optimal thread count (insufficient data).")
        end
        println(io)
        println(io, "## Results Table")
        println(io)
        println(io, "- `time_per_case_s`: wall_time / n_cases (n_cases = workers × batches). Lower = better throughput.")
        println(io, "- `throughput_cases_s`: n_cases / wall_time. Higher = better.")
        println(io, "- `speedup_throughput`: throughput_at_N / throughput_at_1. Ideal = N.")
        println(io, "- `efficiency_pct`: 100 × speedup / N. 100% = perfect linear scaling. Expect a knee at the physical core count; efficiency past it reflects SMT, not real cores.")
        println(io, "- `imbalance_pct`: 100 × (max_case − min_case) / mean_case. High values indicate uneven thread utilisation.")
        println(io, "- `mean_gc_pct`: 100 × Σgc_time / Σcase_time across cases. High values indicate allocation pressure.")
        println(io)
        cols = [
            :julia_threads, :n_cases, :time_per_case_s, :throughput_cases_s,
            :speedup_throughput, :efficiency_pct,
            :wall_time_s, :mean_case_s, :min_case_s, :max_case_s, :imbalance_pct,
            :mean_gc_pct,
            :is_optimal,
        ]
        present = [c for c in cols if c in Symbol.(names(summary_df))]
        println(io, "| ", join(String.(present), " | "), " |")
        println(io, "|", join(fill("---:", length(present)), "|"), "|")
        for row in eachrow(summary_df)
            vals = map(present) do col
                v = row[col]
                v isa Missing && return "n/a"
                v isa Bool    && return v ? "**optimal**" : ""
                v isa Float64 && return isfinite(v) ? string(round(v; digits=4)) : "n/a"
                string(v)
            end
            println(io, "| ", join(vals, " | "), " |")
        end
        println(io)
        println(io, "## Interpretation")
        println(io)
        println(io, "This benchmark measures **outer parallelism** for Monte Carlo workloads using the same")
        println(io, "distributed execution pattern as `aerobraking_perturbation_mc/main.jl`: N distributed")
        println(io, "worker processes are started via `addprocs`, each running one independent single-satellite")
        println(io, "simulation via `remotecall_fetch` with the same work-stealing `@async` loop used by")
        println(io, "`_run_samples_distributed` in study.jl. The metric `time_per_case = wall_time / N`")
        println(io, "reaches its minimum at the optimal worker count.")
        println(io)
        println(io, "Each worker builds its own simulation config locally (including SPICE calls), so there")
        println(io, "is no shared Fortran global state between concurrent cases. Workers run with")
        println(io, "`--threads=1` and `SPACEAGORA_INNER_THREAD_BUDGET=1` to disable inner parallelism.")
        println(io)
        println(io, "Typical behaviour:")
        println(io, "- At low N: throughput scales nearly linearly (efficiency ≈ 100%) because cases run")
        println(io, "  independently with minimal contention.")
        println(io, "- Near the physical core count: memory bandwidth and inter-process communication")
        println(io, "  overhead cause efficiency to fall below 100%.")
        println(io, "- Beyond the physical core count: OS scheduling overhead raises `time_per_case`")
        println(io, "  and the optimal point is clearly visible.")
        println(io)
        println(io, "This benchmark differs from `performance_thread_scaling_1024.jl`, which measures")
        println(io, "**inner** thread scaling (one multi-satellite simulation with more threads). Here all")
        println(io, "simulations are single-satellite and the parallelism is entirely between independent cases.")
        println(io)
        println(io, "## Reproducibility")
        println(io)
        sysimage = _mc_ts_sysimage()
        project  = _mc_ts_project()
        println(io, "```bash")
        for n in thread_counts
            println(io,
                "SPACEAGORA_MC_THREAD_OUTDIR=$(joinpath(dirname(path), "$(n)t")) " *
                "SPACEAGORA_MC_THREAD_NORBITS=$(n_orbits) " *
                "SPACEAGORA_MC_DIST_NWORKERS=$(n) " *
                (sysimage === nothing ? "" : "SPACEAGORA_MC_THREAD_SYSIMAGE=$(sysimage) ") *
                "julia --threads=1 " *
                (sysimage === nothing ? "" : "--sysimage=$(sysimage) ") *
                "--project=$(project) benchmarks/studies/performance_mc_thread_scaling_worker.jl"
            )
        end
        println(io, "```")
        if !isempty(plot_paths)
            println(io)
            println(io, "## Plots")
            println(io)
            for p in plot_paths
                println(io, "- `$(p)`")
            end
        end
    end
end

function _mc_ts_write_optimal_toml(outdir::String, stamp::String, opt_row)
    opt_row === nothing && return
    toml_path = joinpath(outdir, "mc_optimal_threads.toml")
    data = Dict(
        "generated_utc"      => string(now(UTC)),
        "stamp"              => stamp,
        "julia_version"      => string(VERSION),
        "optimal_threads"    => opt_row.julia_threads,
        "time_per_case_s"    => opt_row.time_per_case_s,
        "throughput_cases_s" => opt_row.throughput_cases_s,
        "case_id"            => opt_row.case_id,
        "n_orbits"           => opt_row.n_orbits,
        "host_cpu_threads"   => Sys.CPU_THREADS,
        "host_physical_cores" => _mc_ts_physical_cores()[1],
    )
    open(toml_path, "w") do io
        TOML.print(io, data)
    end
    println("[mc-thread-scaling] optimal thread count written to: $(toml_path)")
end

function main_mc_thread_scaling()
    outdir        = _mc_ts_env("SPACEAGORA_MC_THREAD_OUTDIR", _MC_TS_DEFAULT_OUTDIR)
    n_orbits      = parse(Int, _mc_ts_env("SPACEAGORA_MC_THREAD_NORBITS", "10"))
    smoke         = _mc_ts_bool("SPACEAGORA_MC_THREAD_SMOKE", false)
    thread_counts = smoke ? [1] : _mc_ts_thread_counts()
    project       = _mc_ts_project()
    sysimage      = _mc_ts_sysimage()

    mkpath(outdir)

    physical_cores, physical_exact = _mc_ts_physical_cores()
    println("MC distributed worker scaling study")
    println("Worker counts:    $(join(thread_counts, ", "))")
    println("Orbits per case:  $(n_orbits)")
    println("Batches/worker:   $(_mc_ts_env("SPACEAGORA_MC_THREAD_BATCHES", "1"))")
    println("Physical cores:   $(physical_cores)$(physical_exact ? "" : " (estimate)")")
    println("Host CPU threads: $(Sys.CPU_THREADS)")
    println("Sysimage:         $(sysimage === nothing ? "none" : sysimage)")
    println("Output:           $(outdir)")
    println("Project:          $(project)")
    println("Smoke:            $(smoke)")
    println()

    rows = []
    for n in thread_counts
        result = _mc_ts_run_one(n, n_orbits, project, outdir)
        result !== nothing && push!(rows, result)
    end

    if isempty(rows)
        @warn "[mc-thread-scaling] no results collected."
        return
    end

    summary_df = _mc_ts_build_summary(rows)
    opt_row    = _mc_ts_optimal_row(summary_df)

    stamp       = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    raw_path    = joinpath(outdir, "mc_thread_scaling_raw_$(stamp).csv")
    summ_path   = joinpath(outdir, "mc_thread_scaling_summary_$(stamp).csv")
    report_path = joinpath(outdir, "mc_thread_scaling_report_$(stamp).md")

    plot_paths = _mc_ts_generate_plots(outdir, stamp, summary_df)

    CSV.write(raw_path, summary_df)
    CSV.write(summ_path, summary_df)
    _mc_ts_write_report(report_path, summary_df, thread_counts, n_orbits; plot_paths=plot_paths)
    _mc_ts_write_optimal_toml(outdir, stamp, opt_row)

    println()
    println("MC thread scaling study complete.")
    if opt_row !== nothing
        @printf("Optimal thread count: %d  (time_per_case=%.3f s, throughput=%.3f cases/s)\n",
            opt_row.julia_threads, opt_row.time_per_case_s, opt_row.throughput_cases_s)
    end
    println("Raw data:    $(raw_path)")
    println("Summary:     $(summ_path)")
    println("Report:      $(report_path)")
    for p in plot_paths
        println("Plot:        $(p)")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_mc_thread_scaling()
end
