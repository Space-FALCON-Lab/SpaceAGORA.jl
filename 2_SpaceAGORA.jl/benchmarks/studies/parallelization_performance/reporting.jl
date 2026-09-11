ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots
using Plots.PlotMeasures: mm

function ppc_hardware_snapshot()
    return (
        timestamp_utc=string(now(UTC)),
        machine=gethostname(),
        julia_version=string(VERSION),
        cpu_threads=Sys.CPU_THREADS,
        julia_threads=Threads.nthreads(),
        os=string(Sys.KERNEL),
        arch=string(Sys.ARCH),
        git_commit=ppc_git_commit()
    )
end

function ppc_git_commit()::String
    try
        return strip(read(`git -C $(PPC_REPO_ROOT) rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end

function ppc_write_rows(path::String, rows::Vector{NamedTuple})
    mkpath(dirname(path))
    df = DataFrame(rows)
    CSV.write(path, df)
    return path
end

function ppc_read_optional(path::String)::DataFrame
    isfile(path) || return DataFrame()
    return CSV.read(path, DataFrame)
end

function ppc_summarize(raw::DataFrame, parity::DataFrame)::DataFrame
    nrow(raw) == 0 && return DataFrame()
    grouped = groupby(raw, [:case, :family, :mode, :thread_count, :mc_samples])
    summary = combine(
        grouped,
        :success => (x -> count(identity, skipmissing(x))) => :success_count,
        :success => length => :sample_count,
        :wall_time_s => mean => :wall_time_mean_s,
        :wall_time_s => std => :wall_time_std_s,
        :wall_time_s => (x -> length(x) <= 1 ? 0.0 : 100.0 * std(x) / max(abs(mean(x)), eps(Float64))) => :wall_time_cv_pct,
        :throughput_samples_per_s => mean => :throughput_samples_per_s_mean,
        :execution_scope => (x -> first(skipmissing(x))) => :execution_scope,
        :outer_backend_actual => (x -> first(skipmissing(x))) => :outer_backend_actual
    )
    summary[!, :wall_time_ci95_low_s] = summary.wall_time_mean_s .- 1.96 .* coalesce.(summary.wall_time_std_s, 0.0) ./ sqrt.(summary.sample_count)
    summary[!, :wall_time_ci95_high_s] = summary.wall_time_mean_s .+ 1.96 .* coalesce.(summary.wall_time_std_s, 0.0) ./ sqrt.(summary.sample_count)

    serial = summary[summary.mode .== "serial", [:case, :mc_samples, :wall_time_mean_s]]
    rename!(serial, :wall_time_mean_s => :serial_wall_time_mean_s)
    summary = leftjoin(summary, serial; on=[:case, :mc_samples])
    summary[!, :speedup_vs_serial] = [
        (ismissing(s) || t <= 0.0) ? missing : s / t
        for (s, t) in zip(summary.serial_wall_time_mean_s, summary.wall_time_mean_s)
    ]
    summary[!, :parallel_efficiency] = [
        (v isa Missing || th <= 0) ? missing : v / max(1, th)
        for (v, th) in zip(summary.speedup_vs_serial, summary.thread_count)
    ]

    if nrow(parity) > 0 && "pass" in names(parity)
        parity_ok = combine(groupby(parity, [:case, :mode]), :pass => (x -> all(Bool.(x))) => :parity_pass)
        summary = leftjoin(summary, parity_ok; on=[:case, :mode])
    else
        summary[!, :parity_pass] = fill(missing, nrow(summary))
    end
    summary[!, :headline_eligible] = [
        ok && (p === missing || p === true)
        for (ok, p) in zip(summary.success_count .== summary.sample_count, summary.parity_pass)
    ]
    return summary
end

function ppc_write_report(path::String, cfg::PPCConfig, raw::DataFrame, summary::DataFrame, parity::DataFrame)
    open(path, "w") do io
        println(io, "# SpaceAGORA Parallelization Performance Study")
        println(io)
        println(io, "- Generated UTC: $(now(UTC))")
        println(io, "- Profile: `$(cfg.profile)`")
        println(io, "- Solver mode: `$(cfg.solver_mode)`")
        println(io, "- Repeats: `$(cfg.repeats)`, warmup: `$(cfg.warmup)`")
        println(io, "- Modes: `$(join(cfg.modes, ", "))`")
        println(io, "- Thread counts: `$(join(cfg.threads, ", "))`")
        println(io)
        println(io, "## Outputs")
        println(io)
        println(io, "- Raw rows: `$(nrow(raw))`")
        println(io, "- Summary rows: `$(nrow(summary))`")
        println(io, "- Trajectory parity rows: `$(nrow(parity))`")
        println(io)
        println(io, "## Headline Summary")
        println(io)
        if nrow(summary) == 0
            println(io, "No summary rows were generated.")
        else
            cols = intersect([:case, :mode, :thread_count, :mc_samples, :execution_scope, :outer_backend_actual, :wall_time_mean_s, :throughput_samples_per_s_mean, :speedup_vs_serial, :headline_eligible], Symbol.(names(summary)))
            show(io, MIME"text/plain"(), first(select(summary, cols), min(30, nrow(summary))))
            println(io)
        end
        println(io)
        println(io, "## Trajectory Correctness")
        println(io)
        if nrow(parity) == 0
            println(io, "No trajectory parity rows were generated.")
        else
            cols = intersect([:case, :mode, :pass, :pos_rel_rms, :pos_rel_max, :vel_rel_rms, :vel_rel_max, :event_count_equal, :event_time_abs_max_s], Symbol.(names(parity)))
            show(io, MIME"text/plain"(), select(parity, cols))
            println(io)
        end
    end
    return path
end

function ppc_write_plots(outdir::String, summary::DataFrame, parity::DataFrame)::Vector{String}
    paths = String[]
    try
        ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
        @eval using Plots
        if nrow(summary) > 0 && "speedup_vs_serial" in names(summary)
            df = summary[summary.mode .!= "serial", :]
            if nrow(df) > 0
                p = Plots.scatter(
                    string.(df.case),
                    coalesce.(df.speedup_vs_serial, NaN);
                    group=string.(df.mode),
                    xrotation=35,
                    xlabel="case",
                    ylabel="speedup vs R0",
                    legend=:outertopright,
                    tickfont=Plots.font(8),
                    guidefont=Plots.font(12),
                    legend_font=Plots.font(9),
                    size=(1500, 850),
                    left_margin=18mm,
                    right_margin=34mm,
                    top_margin=10mm,
                    bottom_margin=36mm
                )
                path = joinpath(outdir, "parallelization_speedup_by_case.png")
                Plots.savefig(p, path)
                push!(paths, path)
            end
        end
        if nrow(summary) > 0 && all(name in names(summary) for name in ("speedup_vs_serial", "execution_scope"))
            inner = summary[(summary.mode .!= "serial") .& (summary.execution_scope .== "single_simulation"), :]
            if nrow(inner) > 0
                p = Plots.scatter(
                    string.(inner.case),
                    coalesce.(inner.speedup_vs_serial, NaN);
                    group=string.(inner.mode),
                    xrotation=35,
                    xlabel="case",
                    ylabel="single-simulation speedup vs R0",
                    legend=:outertopright,
                    tickfont=Plots.font(8),
                    guidefont=Plots.font(12),
                    legend_font=Plots.font(9),
                    size=(1500, 850),
                    left_margin=18mm,
                    right_margin=34mm,
                    top_margin=10mm,
                    bottom_margin=36mm
                )
                path = joinpath(outdir, "parallelization_inner_speedup_by_case.png")
                Plots.savefig(p, path)
                push!(paths, path)
            end
        end
        if nrow(summary) > 0 && all(name in names(summary) for name in ("throughput_samples_per_s_mean", "outer_backend_actual"))
            outer = summary[(summary.mode .!= "serial") .& (summary.mc_samples .> 1) .& (summary.outer_backend_actual .!= "serial"), :]
            if nrow(outer) > 0
                p = Plots.scatter(
                    outer.mc_samples,
                    outer.throughput_samples_per_s_mean;
                    group=string.(outer.mode) .* " / " .* string.(outer.case),
                    xscale=:log10,
                    xlabel="Monte Carlo samples",
                    ylabel="throughput (samples/s)",
                    legend=:outertopright,
                    tickfont=Plots.font(8),
                    guidefont=Plots.font(12),
                    legend_font=Plots.font(8),
                    size=(1500, 850),
                    left_margin=18mm,
                    right_margin=48mm,
                    top_margin=10mm,
                    bottom_margin=18mm
                )
                path = joinpath(outdir, "parallelization_outer_throughput_by_samples.png")
                Plots.savefig(p, path)
                push!(paths, path)
            end
        end
        if nrow(parity) > 0 && "pos_rel_max" in names(parity)
            p = Plots.scatter(
                string.(parity.case),
                parity.pos_rel_max;
                group=string.(parity.mode),
                xrotation=45,
                yscale=:log10,
                ylabel="max relative position error",
                legend=:outertopright,
                size=(1200, 700)
            )
            path = joinpath(outdir, "parallelization_parity_position_error.png")
            Plots.savefig(p, path)
            push!(paths, path)
        end
    catch err
        @warn "Plot generation failed" exception=(err, catch_backtrace())
    end
    return paths
end
