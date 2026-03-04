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
    bench_elapsed_s::Float64
    orbit_elapsed_s::Float64
    entry_duration_elapsed_s::Float64
    raw_path::String
    summary_path::String
    orbit_raw_path::String
    orbit_summary_path::String
    entry_duration_raw_path::String
    entry_duration_summary_path::String
    report_path::String
    stage_timing_path::String = ""
    hardware_info_path::String = ""
    split_gate_elapsed_s::Float64 = 0.0
    split_gate_csv_path::Union{Nothing, String} = nothing
    split_gate_report_path::Union{Nothing, String} = nothing
    split_gate_df::Union{Nothing, DataFrame} = nothing
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

    raw_df = DataFrame()
    summary_df = DataFrame()
    orbit_raw_df = DataFrame()
    orbit_summary_df = DataFrame()
    entry_duration_raw_df = DataFrame()
    entry_duration_summary_df = DataFrame()
    density_backend_breakdown_df = DataFrame()
    bench_elapsed_s = 0.0
    split_gate_elapsed_s = 0.0
    split_gate_df = nothing
    split_gate_csv_path = nothing
    split_gate_report_path = nothing
    orbit_elapsed_s = 0.0
    entry_duration_elapsed_s = 0.0
    withenv("SPACEAGORA_PERF_PARALLEL_BACKEND" => backend) do
        bench_started_ns = time_ns()
        raw_df = run_benchmarks(config.profile, cases, planet)
        summary_df = summarize_results(raw_df)
        density_backend_breakdown_df = summarize_density_backend_breakdown(raw_df)
        bench_elapsed_s = (time_ns() - bench_started_ns) / 1e9

        if _split_rollout_enabled()
            split_gate_started_ns = time_ns()
            split_gate_result = evaluate_split_rollout_gate(config.profile, cases, mode_outdir)
            split_gate_elapsed_s = (time_ns() - split_gate_started_ns) / 1e9
            split_gate_df = split_gate_result.df
            split_gate_csv_path = split_gate_result.csv_path
            split_gate_report_path = split_gate_result.report_path
        end

        orbit_started_ns = time_ns()
        orbit_raw_df = run_per_orbit_for_scenarios(config.profile, cases, planet)
        orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)
        orbit_elapsed_s = (time_ns() - orbit_started_ns) / 1e9

        entry_duration_started_ns = time_ns()
        entry_duration_raw_df = run_entry_duration_sweep(config.profile, cases, planet)
        entry_duration_summary_df = summarize_entry_duration_results(entry_duration_raw_df)
        entry_duration_elapsed_s = (time_ns() - entry_duration_started_ns) / 1e9
    end
    elapsed_s = bench_elapsed_s + split_gate_elapsed_s + orbit_elapsed_s + entry_duration_elapsed_s

    raw_path = joinpath(mode_outdir, "runtime_raw_$(config.profile.name)_$(mode)_$(stamp).csv")
    summary_path = joinpath(mode_outdir, "runtime_summary_$(config.profile.name)_$(mode)_$(stamp).csv")
    orbit_raw_path = joinpath(mode_outdir, "runtime_per_orbit_raw_$(config.profile.name)_$(mode)_$(stamp).csv")
    orbit_summary_path = joinpath(mode_outdir, "runtime_per_orbit_summary_$(config.profile.name)_$(mode)_$(stamp).csv")
    entry_duration_raw_path = joinpath(mode_outdir, "runtime_entry_duration_raw_$(config.profile.name)_$(mode)_$(stamp).csv")
    entry_duration_summary_path = joinpath(mode_outdir, "runtime_entry_duration_summary_$(config.profile.name)_$(mode)_$(stamp).csv")
    stage_timing_path = joinpath(mode_outdir, "runtime_stage_timing_$(config.profile.name)_$(mode)_$(stamp).csv")
    hardware_info_path = joinpath(mode_outdir, "runtime_hardware_info_$(config.profile.name)_$(mode)_$(stamp).csv")
    report_path = joinpath(mode_outdir, "runtime_report_$(config.profile.name)_$(mode)_$(stamp).md")
    plot_stamp = "$(mode)_$(stamp)"
    hw = _runtime_hardware_snapshot()
    inner_hint_layer_df = _inner_hint_layer_report_df(config.profile, hw)
    plot_paths = generate_runtime_plots(
        mode_outdir,
        config.profile,
        plot_stamp,
        raw_df,
        summary_df,
        orbit_summary_df,
        entry_duration_summary_df;
        split_gate_df=split_gate_df,
        inner_hint_layer_df=inner_hint_layer_df,
        density_backend_breakdown_df=density_backend_breakdown_df
    )

    stage_names = ["run_benchmarks"]
    stage_elapsed = [bench_elapsed_s]
    if _split_rollout_enabled()
        push!(stage_names, "run_split_rollout_gate")
        push!(stage_elapsed, split_gate_elapsed_s)
    end
    push!(stage_names, "run_per_orbit")
    push!(stage_names, "run_entry_duration_sweep")
    push!(stage_names, "total")
    push!(stage_elapsed, orbit_elapsed_s)
    push!(stage_elapsed, entry_duration_elapsed_s)
    push!(stage_elapsed, elapsed_s)
    stage_timing_df = DataFrame(stage=stage_names, elapsed_s=stage_elapsed)
    hardware_info_df = DataFrame([
        (
            profile=config.profile.name,
            mode=string(mode),
            machine_label=hw.machine_label,
            hardware_class=hw.hardware_class,
            host_name=hw.host_name,
            cpu_name=hw.cpu_name,
            cpu_threads=hw.cpu_threads,
            julia_threads=hw.julia_threads,
            os=hw.os,
            arch=hw.arch
        )
    ])

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    CSV.write(orbit_raw_path, orbit_raw_df)
    CSV.write(orbit_summary_path, orbit_summary_df)
    CSV.write(entry_duration_raw_path, entry_duration_raw_df)
    CSV.write(entry_duration_summary_path, entry_duration_summary_df)
    CSV.write(stage_timing_path, stage_timing_df)
    CSV.write(hardware_info_path, hardware_info_df)
    write_report(
        report_path,
        config.profile,
        raw_df,
        summary_df,
        orbit_summary_df,
        entry_duration_summary_df;
        plot_paths=plot_paths,
        split_gate_df=split_gate_df,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_report_path=split_gate_report_path,
        stage_timing_df=stage_timing_df,
        entry_duration_summary_csv_path=entry_duration_summary_path
    )

    println(
        "[paper-pipeline] mode=$(mode) completed in $(round(elapsed_s; digits=3)) s " *
        "(run_benchmarks=$(round(bench_elapsed_s; digits=3)) s, split_gate=$(round(split_gate_elapsed_s; digits=3)) s, " *
        "mission_time_sweep=$(round(orbit_elapsed_s; digits=3)) s, entry_duration_sweep=$(round(entry_duration_elapsed_s; digits=3)) s, plots=$(length(plot_paths)))"
    )

    return ModeRunArtifacts(
        mode=mode,
        backend=backend,
        elapsed_s=elapsed_s,
        bench_elapsed_s=bench_elapsed_s,
        orbit_elapsed_s=orbit_elapsed_s,
        entry_duration_elapsed_s=entry_duration_elapsed_s,
        raw_path=raw_path,
        summary_path=summary_path,
        orbit_raw_path=orbit_raw_path,
        orbit_summary_path=orbit_summary_path,
        entry_duration_raw_path=entry_duration_raw_path,
        entry_duration_summary_path=entry_duration_summary_path,
        report_path=report_path,
        stage_timing_path=stage_timing_path,
        hardware_info_path=hardware_info_path,
        split_gate_elapsed_s=split_gate_elapsed_s,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_report_path=split_gate_report_path,
        split_gate_df=split_gate_df,
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
        :total_time_mean_s,
        :total_time_ci95_low_s,
        :total_time_ci95_high_s,
        :total_time_cv_pct
    )
    rename!(comparison, :total_time_mean_s => :serial_total_time_mean_s)
    rename!(comparison, :total_time_ci95_low_s => :serial_total_time_ci95_low_s)
    rename!(comparison, :total_time_ci95_high_s => :serial_total_time_ci95_high_s)
    rename!(comparison, :total_time_cv_pct => :serial_total_time_cv_pct)

    for artifact in artifacts
        mode = artifact.mode
        mode == :serial && continue
        mode_df = select(
            artifact.summary_df,
            :scenario,
            :total_time_mean_s,
            :total_time_ci95_low_s,
            :total_time_ci95_high_s,
            :total_time_cv_pct
        )
        mean_col = Symbol("$(mode)_total_time_mean_s")
        ci_low_col = Symbol("$(mode)_total_time_ci95_low_s")
        ci_high_col = Symbol("$(mode)_total_time_ci95_high_s")
        cv_col = Symbol("$(mode)_total_time_cv_pct")
        rename!(mode_df, :total_time_mean_s => mean_col)
        rename!(mode_df, :total_time_ci95_low_s => ci_low_col)
        rename!(mode_df, :total_time_ci95_high_s => ci_high_col)
        rename!(mode_df, :total_time_cv_pct => cv_col)
        comparison = leftjoin(comparison, mode_df, on=:scenario)
    end

    for artifact in artifacts
        mode = artifact.mode
        speedup_col = Symbol("speedup_vs_serial_$(mode)")
        if mode == :serial
            comparison[!, speedup_col] = ones(Float64, nrow(comparison))
        else
            mean_col = Symbol("$(mode)_total_time_mean_s")
            serial_vals = _safe_float.(comparison.serial_total_time_mean_s)
            mode_vals = _safe_float.(comparison[!, mean_col])
            comparison[!, speedup_col] = [
                (!(s isa Missing) && !(m isa Missing) && isfinite(Float64(s)) && Float64(s) > 0.0 && isfinite(Float64(m)) && Float64(m) > 0.0) ? (Float64(s) / Float64(m)) : missing
                for (s, m) in zip(serial_vals, mode_vals)
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
        machine_label = (:machine_label in names(raw_df) && nrow(raw_df) > 0) ? string(raw_df.machine_label[1]) : _machine_label()
        hardware_class = (:hardware_class in names(raw_df) && nrow(raw_df) > 0) ? string(raw_df.hardware_class[1]) : _hardware_class_name()
        failures = count(!, raw_df.solve_success)
        unstable = count(occursin("Unstable"), string.(raw_df.solve_retcode))
        baseline_idx = findfirst(==(PERF_BASELINE_SCENARIO), summary_df.scenario)
        baseline_mean = baseline_idx === nothing ? missing : summary_df.total_time_mean_s[baseline_idx]
        mc = raw_df[(raw_df.category .== "montecarlo") .& (raw_df.solve_success .== true), :]
        mc_mean = nrow(mc) > 0 ? mean(mc.total_time_s) : missing
        mc_p90 = nrow(mc) > 0 ? quantile(mc.total_time_s, 0.9) : missing
        elapsed_total = artifact.elapsed_s
        bench_elapsed = artifact.bench_elapsed_s
        split_elapsed = artifact.split_gate_elapsed_s
        orbit_elapsed = artifact.orbit_elapsed_s
        entry_duration_elapsed = artifact.entry_duration_elapsed_s
        orbit_share = elapsed_total > 0.0 ? (100.0 * orbit_elapsed / elapsed_total) : missing
        entry_duration_share = elapsed_total > 0.0 ? (100.0 * entry_duration_elapsed / elapsed_total) : missing
        split_gate_total = 0
        split_gate_pass = 0
        if !(artifact.split_gate_df === nothing) && (:pass_all in names(artifact.split_gate_df))
            split_gate_total = nrow(artifact.split_gate_df)
            split_gate_pass = count(Bool.(artifact.split_gate_df.pass_all))
        end
        split_gate_pass_rate = split_gate_total > 0 ? (100.0 * split_gate_pass / split_gate_total) : missing

        push!(rows, (
            mode=string(artifact.mode),
            backend=artifact.backend,
            machine_label=machine_label,
            hardware_class=hardware_class,
            elapsed_s=elapsed_total,
            bench_elapsed_s=bench_elapsed,
            split_gate_elapsed_s=split_elapsed,
            orbit_elapsed_s=orbit_elapsed,
            entry_duration_elapsed_s=entry_duration_elapsed,
            orbit_share_pct=orbit_share,
            entry_duration_share_pct=entry_duration_share,
            rows_raw=nrow(raw_df),
            rows_summary=nrow(summary_df),
            failed_rows=failures,
            unstable_rows=unstable,
            baseline_mean_s=baseline_mean,
            montecarlo_mean_s=mc_mean,
            montecarlo_p90_s=mc_p90,
            split_gate_rows=split_gate_total,
            split_gate_pass_rows=split_gate_pass,
            split_gate_pass_rate_pct=split_gate_pass_rate
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

@inline function _pipeline_plot_scenario_limit()::Int
    raw = strip(get(ENV, "SPACEAGORA_PIPELINE_PLOT_SCENARIO_LIMIT", "0"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_PIPELINE_PLOT_SCENARIO_LIMIT must be an integer, got '$raw'"))
    end
    return max(0, parsed)
end

@inline function _pipeline_speedup_series_color(label::String)::String
    token = lowercase(strip(label))
    if token in ("auto", "auto_static", "r2_inner_only", "inner_only")
        return "#2f7fc1"
    elseif token in ("auto_adaptive", "r4_outer_inner_adaptive", "outer_inner_adaptive")
        return "#4f9d69"
    elseif token in ("threads", "threads_static", "r1_outer_only", "r1_a_outer_only", "outer_only")
        return "#d67c1c"
    elseif token in ("r0p_outer_only_process", "r1_b_outer_only_process", "outer_only_process")
        return "#b64a3a"
    elseif token in ("process", "process_static", "r3_outer_inner_static", "outer_inner_static")
        return "#7b63c6"
    end
    # Fallback keeps unknown/new series visible and distinct from the default auto color.
    return "#c44e52"
end

@inline function _pipeline_speedup_series_label(label::String)::String
    token = lowercase(strip(label))
    if token in ("outer_only_process", "r0p_outer_only_process", "r1_b_outer_only_process")
        return "R1_b outer-only process"
    end
    if token in ("outer_only", "r1_outer_only", "r1_a_outer_only")
        return "R1_a outer-only"
    elseif token in ("inner_only", "r2_inner_only")
        return "R2 inner-only"
    elseif token in ("outer_inner_static", "r3_outer_inner_static")
        return "R3 outer+inner static"
    elseif token in ("outer_inner_adaptive", "r4_outer_inner_adaptive")
        return "R4 outer+inner adaptive"
    elseif token in ("auto", "auto_static")
        return "Auto static"
    elseif token == "auto_adaptive"
        return "Auto adaptive"
    elseif token in ("threads", "threads_static")
        return "Threads"
    elseif token in ("process", "process_static")
        return "Process"
    end
    return label
end

@inline function _pipeline_speedup_palette(series_labels::Vector{String})::Vector{String}
    return [_pipeline_speedup_series_color(label) for label in series_labels]
end

@inline function _pipeline_series_visible(vals::AbstractVector{<:Real}, eps::Float64)::Bool
    return any(isfinite(v) && abs(v) > eps for v in vals)
end

@inline function _pipeline_series_visible(vals::AbstractVector{<:Real})::Bool
    return _pipeline_series_visible(vals, 1e-9)
end

@inline function _pipeline_speedup_series_rank(token::String)::Int
    t = lowercase(strip(token))
    if t in ("serial", "r0_true_serial")
        return 0
    elseif t in ("r1_b_outer_only_process", "outer_only_process", "r0p_outer_only_process")
        return 1
    elseif t in ("r1_a_outer_only", "outer_only", "r1_outer_only", "threads", "threads_static")
        return 2
    elseif t in ("r2_inner_only", "inner_only", "auto", "auto_static")
        return 3
    elseif t in ("r3_outer_inner_static", "outer_inner_static", "process", "process_static")
        return 4
    elseif t in ("r4_outer_inner_adaptive", "outer_inner_adaptive", "auto_adaptive")
        return 5
    end
    return typemax(Int)
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

function _pipeline_top_scenarios(comparison_df::DataFrame; limit::Int=0)::DataFrame
    row_limit = limit <= 0 ? nrow(comparison_df) : limit
    if !_has_column(comparison_df, :serial_total_time_mean_s)
        return first(comparison_df, min(row_limit, nrow(comparison_df)))
    end
    if limit <= 0
        return copy(comparison_df)
    end
    df = copy(comparison_df)
    df[!, :_serial_sort_key] = [ismissing(v) ? -Inf : Float64(v) for v in df.serial_total_time_mean_s]
    sort!(df, :_serial_sort_key, rev=true)
    select!(df, Not(:_serial_sort_key))
    return first(df, min(row_limit, nrow(df)))
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

    # 2) Stage-split elapsed by mode.
    if nrow(overview_df) > 0 && _has_column(overview_df, :bench_elapsed_s) && _has_column(overview_df, :orbit_elapsed_s)
        bench = Float64.(overview_df.bench_elapsed_s)
        split_gate = _has_column(overview_df, :split_gate_elapsed_s) ? Float64.(overview_df.split_gate_elapsed_s) : zeros(length(bench))
        orbit = Float64.(overview_df.orbit_elapsed_s)
        entry_duration = _has_column(overview_df, :entry_duration_elapsed_s) ? Float64.(overview_df.entry_duration_elapsed_s) : zeros(length(bench))
        stage_values = [bench, split_gate, orbit, entry_duration]
        stage_labels = ["run_benchmarks", "split_rollout_gate", "per_orbit", "entry_duration_sweep"]
        stage_colors = ["#4b86b4", "#6b5ca5", "#d97706", "#2f855a"]
        stage_maxima = [maximum([abs(v) for v in vals if isfinite(v)]; init=0.0) for vals in stage_values]
        global_stage_max = maximum(stage_maxima; init=0.0)
        stage_eps = max(1e-9, 1e-4 * global_stage_max)
        keep = [_pipeline_series_visible(vals, stage_eps) for vals in stage_values]
        if any(keep)
            selected_values = stage_values[keep]
            vals = hcat(selected_values...)
            labels = stage_labels[keep]
            colors = stage_colors[keep]
            label_value = length(labels) == 1 ? labels[1] : labels
            color_value = length(colors) == 1 ? colors[1] : colors
            plt = Plots.bar(
                mode_labels,
                vals;
                label=label_value,
                color=color_value,
                bar_position=:stack,
                title="Pipeline Stage Wall Time by Mode",
                xlabel="Mode",
                ylabel="Elapsed (s)",
                _plot_margins(size=(1900, 1000), bottom_mm=20, right_mm=50, legend=:outertopright)...
            )
            _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_mode_stage_elapsed", profile, stamp)
        end
    end

    # 3) Baseline scenario mean runtime by mode.
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
                title="Baseline Runtime by Mode ($(PERF_BASELINE_SCENARIO))",
                xlabel="Mode",
                ylabel="Mean total time (s)",
                _plot_margins(size=(1700, 950), bottom_mm=18, right_mm=18)...
            )
            _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_mode_baseline", profile, stamp)
        end
    end

    # 4) Monte Carlo p90 runtime by mode.
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

    scenario_limit = _pipeline_plot_scenario_limit()
    top_df = _pipeline_top_scenarios(comparison_df; limit=scenario_limit)
    scenario_scope = scenario_limit <= 0 ? "All Scenarios" : "Top Scenarios by Serial Cost"
    scenario_labels = [_plot_axis_label(String(s)) for s in top_df.scenario]
    x = collect(1:length(scenario_labels))

    # 5) Speedup vs serial by scenario (for each non-serial mode).
    speedup_cols = [name for name in names(top_df) if startswith(String(name), "speedup_vs_serial_") && String(name) != "speedup_vs_serial_serial"]
    if !isempty(speedup_cols) && !isempty(scenario_labels)
        sort!(speedup_cols; by=c -> (_pipeline_speedup_series_rank(replace(String(c), "speedup_vs_serial_" => "")), String(c)))
        values = Matrix{Float64}(undef, nrow(top_df), length(speedup_cols))
        for (j, col) in enumerate(speedup_cols)
            values[:, j] = [
                (v isa Missing || !isfinite(Float64(v))) ? NaN : Float64(v)
                for v in top_df[!, col]
            ]
        end
        finite_series = [any(isfinite.(values[:, j])) for j in 1:size(values, 2)]
        if !any(finite_series)
            finite_series .= true
        end
        values = values[:, finite_series]
        kept_cols = speedup_cols[finite_series]
        series_tokens = [replace(String(c), "speedup_vs_serial_" => "") for c in kept_cols]
        series_labels = [_pipeline_speedup_series_label(token) for token in series_tokens]
        series_palette = _pipeline_speedup_palette(series_tokens)
        label_value = length(series_labels) == 1 ? series_labels[1] : series_labels
        color_value = length(series_palette) == 1 ? series_palette[1] : series_palette
        plt = Plots.bar(
            x,
            values;
            label=label_value,
            color=color_value,
            title="Speedup vs Serial ($scenario_scope)",
            xlabel="Scenario",
            ylabel="Speedup (x)",
            xticks=(x, scenario_labels),
            _plot_margins(size=(2400, 1200), bottom_mm=92, right_mm=62, legend=:outertopright)...
        )
        Plots.hline!(plt, [1.0]; color=:black, linestyle=:dash, label="Serial = 1x")
        _save_pipeline_plot!(plot_paths, plt, outdir, "paper_plot_speedup_vs_serial", profile, stamp)
    end

    # 6) Scenario runtime heatmap across modes.
    runtime_cols = [name for name in names(top_df) if endswith(String(name), "_total_time_mean_s")]
    if !isempty(runtime_cols) && !isempty(scenario_labels)
        mode_order = [replace(String(c), "_total_time_mean_s" => "") for c in runtime_cols]
        scenario_labels_matrix = [_plot_wrapped_label(String(s); width=22, max_lines=3) for s in top_df.scenario]
        matrix_height = max(1300, 420 + 72 * length(scenario_labels_matrix))
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
            yticks=(1:length(scenario_labels), scenario_labels_matrix),
            colorbar=true,
            colorbar_title="mean total (s)",
            xlabel="Mode",
            ylabel="Scenario",
            title="Scenario Runtime Matrix ($scenario_scope)",
            color=Plots.cgrad([:lightsteelblue1, :mediumpurple3, "#3a0a2a"]),
            _plot_margins(size=(2400, matrix_height), left_mm=120, right_mm=32, bottom_mm=18)...
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
    unique_hardware = (:hardware_class in names(overview_df)) ? unique(String.(overview_df.hardware_class)) : String[]
    unique_machines = (:machine_label in names(overview_df)) ? unique(String.(overview_df.machine_label)) : String[]
    total_split_rows = 0
    total_split_pass = 0
    for artifact in artifacts
        if !(artifact.split_gate_df === nothing) && (:pass_all in names(artifact.split_gate_df))
            total_split_rows += nrow(artifact.split_gate_df)
            total_split_pass += count(Bool.(artifact.split_gate_df.pass_all))
        end
    end

    open(path, "w") do io
        println(io, "# SpaceAGORA Paper Runtime Pipeline")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Threads in this Julia process: `$nthreads`")
        println(io, "- Detected CPU threads: `$cpu_threads`")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Modes executed: `$mode_list`")
        if !isempty(unique_hardware)
            println(io, "- Hardware classes observed: `$(join(unique_hardware, ", "))`")
        end
        if !isempty(unique_machines)
            println(io, "- Machine labels observed: `$(join(unique_machines, ", "))`")
        end
        println(io)
        println(io, "## Claim Scope")
        println(io)
        if config.profile.name == "full"
            println(io, "- `full` profile outputs are intended for main paper claims.")
        else
            println(io, "- `quick` profile outputs are development/CI/regression evidence only; not main paper-claim evidence.")
        end
        println(io, "- Mission-time sweep excludes `entry` scenarios by design (baseline-orbit multipliers are not entry-comparable).")
        println(io, "- Entry behavior is represented by the dedicated entry-duration sweep artifacts per mode.")
        println(io)

        println(io, "## Mode Overview")
        println(io)
        _write_markdown_table(io, overview_df)
        println(io)

        println(io, "## Scenario Comparison (Means + CI)")
        println(io)
        _write_markdown_table(io, comparison_df)
        println(io)

        println(io, "## Verification Snapshot")
        println(io)
        if total_split_rows > 0
            pass_rate = 100.0 * total_split_pass / total_split_rows
            println(io, "- Split rollout gate rows: `$(total_split_rows)`; pass rows: `$(total_split_pass)` (`$(round(pass_rate; digits=2))%`).")
        else
            println(io, "- Split rollout gate rows: none (gate disabled or not configured).")
        end
        if :failed_rows in names(overview_df)
            total_failed = sum(Int.(overview_df.failed_rows))
            total_rows = sum(Int.(overview_df.rows_raw))
            println(io, "- Solver-success samples across modes: `$(total_rows - total_failed)/$(total_rows)`.")
        end
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
                baseline_idx = findfirst(==(PERF_BASELINE_SCENARIO), serial_summary.scenario)
                if baseline_idx === nothing
                    println(io, "Baseline scenario `$(PERF_BASELINE_SCENARIO)` not present.")
                else
                    baseline = serial_summary.total_time_mean_s[baseline_idx]
                    println(io, "- Baseline (`$(PERF_BASELINE_SCENARIO)`) mean total time: `$(round(baseline; digits=6)) s`.")
                    for scenario in (
                        "single_orientation_aero",
                        "single_entry_earth_shallow",
                        "single_entry_earth_nominal",
                        "single_entry_earth_steep",
                        "single_nbody_sun_moon",
                        "single_harmonics_l20",
                        "single_harmonics_l50",
                        "single_baseline_long_mission",
                        "multi_4_gravity",
                        "multi_8_gravity",
                        "multi_16_gravity",
                        "multi_32_gravity",
                        "multi_64_gravity"
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
            println(io, "- run_benchmarks elapsed: `$(round(artifact.bench_elapsed_s; digits=3)) s`")
            println(io, "- split rollout gate elapsed: `$(round(artifact.split_gate_elapsed_s; digits=3)) s`")
            println(io, "- mission-time-sweep elapsed: `$(round(artifact.orbit_elapsed_s; digits=3)) s`")
            println(io, "- entry-duration-sweep elapsed: `$(round(artifact.entry_duration_elapsed_s; digits=3)) s`")
            println(io, "- Raw: `$(artifact.raw_path)`")
            println(io, "- Summary: `$(artifact.summary_path)`")
            println(io, "- Mission-time-sweep raw: `$(artifact.orbit_raw_path)`")
            println(io, "- Mission-time-sweep summary: `$(artifact.orbit_summary_path)`")
            println(io, "- Entry-duration-sweep raw: `$(artifact.entry_duration_raw_path)`")
            println(io, "- Entry-duration-sweep summary: `$(artifact.entry_duration_summary_path)`")
            if !isempty(artifact.stage_timing_path)
                println(io, "- Stage timing: `$(artifact.stage_timing_path)`")
            end
            if !isempty(artifact.hardware_info_path)
                println(io, "- Hardware info: `$(artifact.hardware_info_path)`")
            end
            if !(artifact.split_gate_csv_path === nothing)
                println(io, "- Split rollout gate CSV: `$(artifact.split_gate_csv_path)`")
            end
            if !(artifact.split_gate_report_path === nothing)
                println(io, "- Split rollout gate report: `$(artifact.split_gate_report_path)`")
            end
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
