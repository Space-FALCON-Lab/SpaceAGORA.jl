const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_PIPELINE_OUTPUT_DIR = joinpath(REPO_ROOT, "output", "performance", "paper_pipeline")

using CSV
using DataFrames
using Dates
using Statistics

include(joinpath(@__DIR__, "performance_runtime_analysis.jl"))

Base.@kwdef struct PipelineConfig
    profile::ProfileSpec
    modes::Vector{Symbol}
    outdir::String
end

Base.@kwdef struct ModeRunArtifacts
    mode::Symbol
    backend::String
    elapsed_s::Float64
    raw_path::String
    summary_path::String
    orbit_raw_path::String
    orbit_summary_path::String
    report_path::String
    raw_df::DataFrame
    summary_df::DataFrame
    orbit_summary_df::DataFrame
end

@inline function _pipeline_backend_for_mode(mode::Symbol)::String
    if mode == :serial
        return "none"
    elseif mode == :auto
        return "auto"
    elseif mode == :threads
        return "threads"
    elseif mode == :process
        return "process"
    end
    throw(ArgumentError("Unsupported mode '$mode'. Valid: serial, auto, threads, process."))
end

@inline function _parse_bool_token(raw::String)::Bool
    token = lowercase(strip(raw))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean token '$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

function _parse_modes(raw::AbstractString)::Vector{Symbol}
    isempty(strip(raw)) && throw(ArgumentError("Mode list cannot be empty."))
    seen = Set{Symbol}()
    modes = Symbol[]
    for token in split(raw, ",")
        mode = Symbol(lowercase(strip(token)))
        mode in (:serial, :auto, :threads, :process) || throw(ArgumentError("Unsupported mode '$mode'. Valid: serial, auto, threads, process."))
        if !(mode in seen)
            push!(modes, mode)
            push!(seen, mode)
        end
    end
    isempty(modes) && throw(ArgumentError("No valid modes were parsed from '$raw'."))
    return modes
end

function parse_pipeline_cli()::PipelineConfig
    profile_name = lowercase(strip(get(ENV, "SPACEAGORA_PAPER_PROFILE", "full")))
    outdir = get(ENV, "SPACEAGORA_PAPER_OUTDIR", DEFAULT_PIPELINE_OUTPUT_DIR)
    modes_raw = get(ENV, "SPACEAGORA_PAPER_MODES", "serial,auto")

    for arg in ARGS
        if arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(strip(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--modes=")
            modes_raw = split(arg, "=", limit=2)[2]
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..., --modes=serial,auto,threads,process."
            ))
        end
    end

    return PipelineConfig(
        profile=_profile_from_name(profile_name),
        modes=_parse_modes(modes_raw),
        outdir=abspath(outdir)
    )
end

function run_mode(
    mode::Symbol,
    config::PipelineConfig,
    planet::Earth,
    cases::Vector{BenchmarkCase}
)::ModeRunArtifacts
    backend = _pipeline_backend_for_mode(mode)
    mode_outdir = joinpath(config.outdir, string(mode))
    mkpath(mode_outdir)
    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")

    println()
    println("[paper-pipeline] mode=$(mode) backend=$(backend)")
    mode_start = time_ns()

    raw_df = DataFrame()
    summary_df = DataFrame()
    orbit_raw_df = DataFrame()
    orbit_summary_df = DataFrame()
    withenv("SPACEAGORA_PERF_PARALLEL_BACKEND" => backend) do
        raw_df = run_benchmarks(config.profile, cases, planet)
        summary_df = summarize_results(raw_df)
        orbit_raw_df = run_per_orbit_for_scenarios(config.profile, cases, planet)
        orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)
    end
    elapsed_s = (time_ns() - mode_start) / 1e9

    raw_path = joinpath(mode_outdir, "runtime_raw_$(config.profile.name)_$(mode)_$(stamp).csv")
    summary_path = joinpath(mode_outdir, "runtime_summary_$(config.profile.name)_$(mode)_$(stamp).csv")
    orbit_raw_path = joinpath(mode_outdir, "runtime_per_orbit_raw_$(config.profile.name)_$(mode)_$(stamp).csv")
    orbit_summary_path = joinpath(mode_outdir, "runtime_per_orbit_summary_$(config.profile.name)_$(mode)_$(stamp).csv")
    report_path = joinpath(mode_outdir, "runtime_report_$(config.profile.name)_$(mode)_$(stamp).md")
    plot_stamp = "$(mode)_$(stamp)"
    plot_paths = generate_runtime_plots(mode_outdir, config.profile, plot_stamp, raw_df, summary_df, orbit_summary_df)

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    CSV.write(orbit_raw_path, orbit_raw_df)
    CSV.write(orbit_summary_path, orbit_summary_df)
    write_report(report_path, config.profile, raw_df, summary_df, orbit_summary_df; plot_paths=plot_paths)

    println("[paper-pipeline] mode=$(mode) completed in $(round(elapsed_s; digits=3)) s (plots=$(length(plot_paths)))")

    return ModeRunArtifacts(
        mode=mode,
        backend=backend,
        elapsed_s=elapsed_s,
        raw_path=raw_path,
        summary_path=summary_path,
        orbit_raw_path=orbit_raw_path,
        orbit_summary_path=orbit_summary_path,
        report_path=report_path,
        raw_df=raw_df,
        summary_df=summary_df,
        orbit_summary_df=orbit_summary_df
    )
end

@inline function _safe_float(x)
    return x isa Missing ? missing : Float64(x)
end

@inline function _fmt_md(v)
    if v isa Missing
        return "n/a"
    elseif v isa Float64
        return isfinite(v) ? string(round(v; digits=6)) : "n/a"
    elseif v isa Bool
        return v ? "true" : "false"
    else
        return string(v)
    end
end

function _write_markdown_table(io, df::DataFrame)
    names_vec = String.(names(df))
    println(io, "| ", join(names_vec, " | "), " |")
    println(io, "|", join(fill("---", length(names_vec)), "|"), "|")
    for row in eachrow(df)
        values = [_fmt_md(row[name]) for name in names(df)]
        println(io, "| ", join(values, " | "), " |")
    end
end

function build_comparison_table(artifacts::Vector{ModeRunArtifacts})::DataFrame
    mode_order = [a.mode for a in artifacts]
    has_serial = :serial in mode_order
    serial_artifact = has_serial ? artifacts[findfirst(==( :serial), mode_order)] : artifacts[1]

    comparison = select(
        serial_artifact.summary_df,
        :scenario,
        :category,
        :samples_success,
        :samples_total,
        :total_time_mean_s
    )
    rename!(comparison, :total_time_mean_s => :serial_total_time_mean_s)

    for artifact in artifacts
        mode = artifact.mode
        mode == :serial && continue
        mode_df = select(artifact.summary_df, :scenario, :total_time_mean_s)
        mean_col = Symbol("$(mode)_total_time_mean_s")
        rename!(mode_df, :total_time_mean_s => mean_col)
        comparison = leftjoin(comparison, mode_df, on=:scenario)
    end

    for artifact in artifacts
        mode = artifact.mode
        speedup_col = Symbol("speedup_vs_serial_$(mode)")
        if mode == :serial
            comparison[!, speedup_col] = ones(Float64, nrow(comparison))
        else
            mean_col = Symbol("$(mode)_total_time_mean_s")
            comparison[!, speedup_col] = [
                (isfinite(s) && s > 0.0 && isfinite(m) && m > 0.0) ? s / m : missing
                for (s, m) in zip(Float64.(comparison.serial_total_time_mean_s), _safe_float.(comparison[!, mean_col]))
            ]
        end
    end

    sort!(comparison, :serial_total_time_mean_s, rev=true)
    return comparison
end

function build_mode_overview(artifacts::Vector{ModeRunArtifacts})::DataFrame
    rows = NamedTuple[]
    for artifact in artifacts
        raw_df = artifact.raw_df
        summary_df = artifact.summary_df
        failures = count(!, raw_df.solve_success)
        unstable = count(occursin("Unstable"), string.(raw_df.solve_retcode))
        baseline_idx = findfirst(==("single_baseline_gravity"), summary_df.scenario)
        baseline_mean = baseline_idx === nothing ? missing : summary_df.total_time_mean_s[baseline_idx]
        mc = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
        mc_mean = nrow(mc) > 0 ? mean(mc.total_time_s) : missing
        mc_p90 = nrow(mc) > 0 ? quantile(mc.total_time_s, 0.9) : missing

        push!(rows, (
            mode=string(artifact.mode),
            backend=artifact.backend,
            elapsed_s=artifact.elapsed_s,
            rows_raw=nrow(raw_df),
            rows_summary=nrow(summary_df),
            failed_rows=failures,
            unstable_rows=unstable,
            baseline_mean_s=baseline_mean,
            montecarlo_mean_s=mc_mean,
            montecarlo_p90_s=mc_p90
        ))
    end
    return DataFrame(rows)
end

@inline function _pipeline_plot_ready()::Bool
    return myid() == 1 && isdefined(Main, :Plots)
end

@inline function _has_column(df::DataFrame, col::Symbol)::Bool
    return col in Symbol.(names(df))
end

function _save_pipeline_plot!(
    paths::Vector{String},
    plt,
    outdir::String,
    basename::String,
    profile::ProfileSpec,
    stamp::String
)
    path = joinpath(outdir, "$(basename)_$(profile.name)_$(stamp).png")
    try
        Plots.savefig(plt, path)
        push!(paths, path)
    catch err
        @warn "[paper-pipeline] failed to save plot $basename: $(_perf_error_text(err))"
    end
    return nothing
end

function _pipeline_top_scenarios(comparison_df::DataFrame; limit::Int=12)::DataFrame
    if !_has_column(comparison_df, :serial_total_time_mean_s)
        return first(comparison_df, min(limit, nrow(comparison_df)))
    end
    df = comparison_df[.!ismissing.(comparison_df.serial_total_time_mean_s), :]
    nrow(df) == 0 && return first(comparison_df, min(limit, nrow(comparison_df)))
    sort!(df, :serial_total_time_mean_s, rev=true)
    return first(df, min(limit, nrow(df)))
end

function generate_pipeline_comparison_plots(
    outdir::String,
    profile::ProfileSpec,
    stamp::String,
    overview_df::DataFrame,
    comparison_df::DataFrame
)::Vector{String}
    !_pipeline_plot_ready() && return String[]
    _ensure_runtime_plot_theme!()
    mkpath(outdir)

    plot_paths = String[]

    mode_labels = String.(overview_df.mode)

    # 1) Wall-clock elapsed by mode.
    if nrow(overview_df) > 0 && _has_column(overview_df, :elapsed_s)
        elapsed = Float64.(overview_df.elapsed_s)
        plt = Plots.bar(
            mode_labels,
            elapsed;
            color="#5b8fb9",
            legend=false,
            title="Pipeline Wall Time by Mode",
            xlabel="Mode",
            ylabel="Elapsed (s)",
            _plot_margins(size=(1700, 950), bottom_mm=18, right_mm=18)...
        )
        _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_mode_elapsed", profile, stamp)
    end

    # 2) Baseline scenario mean runtime by mode.
    if nrow(overview_df) > 0 && _has_column(overview_df, :baseline_mean_s)
        baseline_vals = Float64[]
        labels = String[]
        for row in eachrow(overview_df)
            val = row.baseline_mean_s
            if !(val isa Missing)
                push!(labels, string(row.mode))
                push!(baseline_vals, Float64(val))
            end
        end
        if !isempty(baseline_vals)
            plt = Plots.bar(
                labels,
                baseline_vals;
                color="#d67c1c",
                legend=false,
                title="Baseline Runtime by Mode (single_baseline_gravity)",
                xlabel="Mode",
                ylabel="Mean total time (s)",
                _plot_margins(size=(1700, 950), bottom_mm=18, right_mm=18)...
            )
            _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_mode_baseline", profile, stamp)
        end
    end

    # 3) Monte Carlo p90 runtime by mode.
    if nrow(overview_df) > 0 && _has_column(overview_df, :montecarlo_p90_s)
        mc_vals = Float64[]
        labels = String[]
        for row in eachrow(overview_df)
            val = row.montecarlo_p90_s
            if !(val isa Missing)
                push!(labels, string(row.mode))
                push!(mc_vals, Float64(val))
            end
        end
        if !isempty(mc_vals)
            plt = Plots.bar(
                labels,
                mc_vals;
                color="#7b63c6",
                legend=false,
                title="Monte Carlo p90 Runtime by Mode",
                xlabel="Mode",
                ylabel="p90 total time (s)",
                _plot_margins(size=(1700, 950), bottom_mm=18, right_mm=18)...
            )
            _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_mode_montecarlo_p90", profile, stamp)
        end
    end

    top_df = _pipeline_top_scenarios(comparison_df; limit=12)
    scenario_labels = [_plot_axis_label(String(s)) for s in top_df.scenario]
    x = collect(1:length(scenario_labels))

    # 4) Speedup vs serial by scenario (for each non-serial mode).
    speedup_cols = [name for name in names(top_df) if startswith(String(name), "speedup_vs_serial_") && String(name) != "speedup_vs_serial_serial"]
    if !isempty(speedup_cols) && !isempty(scenario_labels)
        values = Matrix{Float64}(undef, nrow(top_df), length(speedup_cols))
        for (j, col) in enumerate(speedup_cols)
            values[:, j] = [
                (v isa Missing || !isfinite(Float64(v))) ? NaN : Float64(v)
                for v in top_df[!, col]
            ]
        end
        series_labels = [replace(String(c), "speedup_vs_serial_" => "") for c in speedup_cols]
        label_value = length(series_labels) == 1 ? series_labels[1] : series_labels
        plt = Plots.bar(
            x,
            values;
            label=label_value,
            color="#4f9d69",
            title="Speedup vs Serial (Top Scenarios by Serial Cost)",
            xlabel="Scenario",
            ylabel="Speedup (x)",
            xticks=(x, scenario_labels),
            _plot_margins(size=(2400, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        Plots.hline!(plt, [1.0]; color=:black, linestyle=:dash, label="Serial = 1x")
        _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_speedup_vs_serial", profile, stamp)
    end

    # 5) Scenario runtime heatmap across modes.
    runtime_cols = [name for name in names(top_df) if endswith(String(name), "_total_time_mean_s")]
    if !isempty(runtime_cols) && !isempty(scenario_labels)
        mode_order = [replace(String(c), "_total_time_mean_s" => "") for c in runtime_cols]
        z = fill(NaN, length(scenario_labels), length(runtime_cols))
        for i in 1:nrow(top_df)
            for (j, col) in enumerate(runtime_cols)
                v = top_df[i, col]
                if !(v isa Missing)
                    z[i, j] = Float64(v)
                end
            end
        end
        plt = Plots.heatmap(
            1:length(runtime_cols),
            1:length(scenario_labels),
            z;
            xticks=(1:length(runtime_cols), mode_order),
            yticks=(1:length(scenario_labels), scenario_labels),
            colorbar=true,
            colorbar_title="mean total (s)",
            xlabel="Mode",
            ylabel="Scenario",
            title="Scenario Runtime Matrix (Top Serial-Cost Scenarios)",
            color=Plots.cgrad([:lightsteelblue1, :mediumpurple3, "#3a0a2a"]),
            _plot_margins(size=(2200, 1300), left_mm=72, right_mm=28, bottom_mm=18)...
        )
        _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_runtime_matrix", profile, stamp)
    end

    return plot_paths
end

function write_pipeline_report(
    path::String,
    config::PipelineConfig,
    overview_df::DataFrame,
    comparison_df::DataFrame,
    artifacts::Vector{ModeRunArtifacts};
    comparison_plot_paths::Vector{String}=String[]
)
    generated = string(now(UTC))
    nthreads = Threads.nthreads()
    cpu_threads = Sys.CPU_THREADS
    mode_list = join(string.(config.modes), ", ")

    open(path, "w") do io
        println(io, "# SpaceAGORA Paper Runtime Pipeline")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Threads in this Julia process: `$nthreads`")
        println(io, "- Detected CPU threads: `$cpu_threads`")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Modes executed: `$mode_list`")
        println(io)

        println(io, "## Mode Overview")
        println(io)
        _write_markdown_table(io, overview_df)
        println(io)

        println(io, "## Scenario Comparison (Means)")
        println(io)
        _write_markdown_table(io, comparison_df)
        println(io)

        println(io, "## Comparison Plots")
        println(io)
        println(io, "- Plot artifacts generated: `$(length(comparison_plot_paths))`")
        if !isempty(comparison_plot_paths)
            for plot_path in comparison_plot_paths
                println(io, "- `$(plot_path)`")
            end
        end
        println(io)

        if :serial in config.modes
            println(io, "## Serial Computational Cost")
            println(io)
            serial_idx = findfirst(==("serial"), overview_df.mode)
            if serial_idx === nothing
                println(io, "Serial mode was requested but not found in overview.")
            else
                serial_summary = artifacts[findfirst(a -> a.mode == :serial, artifacts)].summary_df
                baseline_idx = findfirst(==("single_baseline_gravity"), serial_summary.scenario)
                if baseline_idx === nothing
                    println(io, "Baseline scenario `single_baseline_gravity` not present.")
                else
                    baseline = serial_summary.total_time_mean_s[baseline_idx]
                    println(io, "- Baseline (`single_baseline_gravity`) mean total time: `$(round(baseline; digits=6)) s`.")
                    for scenario in (
                        "single_orientation_aero",
                        "single_nbody_sun_moon",
                        "single_harmonics_l20",
                        "single_harmonics_l50",
                        "single_thruster_control",
                        "single_baseline_long_mission",
                        "multi_2_gravity",
                        "multi_4_gravity",
                        "multi_8_gravity",
                        "single_j2"
                    )
                        idx = findfirst(==(scenario), serial_summary.scenario)
                        if idx !== nothing
                            ratio = serial_summary.total_time_mean_s[idx] / baseline
                            println(io, "- `$(scenario)`: `$(round(ratio; digits=3))x` baseline.")
                        end
                    end
                end
            end
            println(io)
        end

        println(io, "## Mode Artifacts")
        println(io)
        for artifact in artifacts
            println(io, "### Mode `$(artifact.mode)` (backend=`$(artifact.backend)`)")
            println(io, "- Wall elapsed: `$(round(artifact.elapsed_s; digits=3)) s`")
            println(io, "- Raw: `$(artifact.raw_path)`")
            println(io, "- Summary: `$(artifact.summary_path)`")
            println(io, "- Per-orbit raw: `$(artifact.orbit_raw_path)`")
            println(io, "- Per-orbit summary: `$(artifact.orbit_summary_path)`")
            println(io, "- Detailed report: `$(artifact.report_path)`")
            println(io)
        end

        println(io, "## Reproducibility")
        println(io)
        println(io, "```bash")
        println(
            io,
            "julia --project=. test/performance_paper_pipeline.jl --profile=$(config.profile.name) --modes=$(join(string.(config.modes), ',')) --outdir=$(config.outdir)"
        )
        println(io, "```")
    end
end

function main()
    config = parse_pipeline_cli()
    mkpath(config.outdir)

    println("Paper pipeline profile: $(config.profile.name)")
    println("Modes: $(join(string.(config.modes), ", "))")
    println("Output directory: $(config.outdir)")

    planet = Earth("", SPICE_PATH)
    cases = build_cases(config.profile, planet)

    artifacts = ModeRunArtifacts[]
    for mode in config.modes
        push!(artifacts, run_mode(mode, config, planet, cases))
    end

    comparison_df = build_comparison_table(artifacts)
    overview_df = build_mode_overview(artifacts)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    comparison_path = joinpath(config.outdir, "paper_compare_summary_$(config.profile.name)_$(stamp).csv")
    overview_path = joinpath(config.outdir, "paper_mode_overview_$(config.profile.name)_$(stamp).csv")
    report_path = joinpath(config.outdir, "paper_pipeline_report_$(config.profile.name)_$(stamp).md")
    comparison_plot_paths = generate_pipeline_comparison_plots(config.outdir, config.profile, stamp, overview_df, comparison_df)

    CSV.write(comparison_path, comparison_df)
    CSV.write(overview_path, overview_df)
    write_pipeline_report(
        report_path,
        config,
        overview_df,
        comparison_df,
        artifacts;
        comparison_plot_paths=comparison_plot_paths
    )

    println()
    println("Pipeline complete.")
    println("Comparison summary: $comparison_path")
    println("Mode overview: $overview_path")
    println("Comparison plots: $(length(comparison_plot_paths))")
    println("Pipeline report: $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
