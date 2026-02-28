const _WRAPPER_DIR = @__DIR__
include(joinpath(_WRAPPER_DIR, "performance_paper_pipeline.jl"))

const STATIC_VS_PARALLEL_PERF_ROOT = joinpath(REPO_ROOT, "output", "performance")
const STATIC_VS_PARALLEL_DEFAULT_OUTDIR = joinpath(STATIC_VS_PARALLEL_PERF_ROOT, "static_vs_parallel")
const STATIC_VS_PARALLEL_RUNTIME_SCRIPT = joinpath(_WRAPPER_DIR, "performance_runtime_analysis.jl")
const STATIC_VS_PARALLEL_PROJECT = joinpath(REPO_ROOT, ".AGORA")

Base.@kwdef struct StaticVsParallelConfig
    profile::ProfileSpec
    outdir::String
    clean::Bool
    include_process::Bool
    process_workers::Union{Nothing, Int}
end

Base.@kwdef struct ArmSpec
    mode::Symbol
    label::String
    backend::String
    adaptive::Bool
    force_inner_off::Bool
    process_workers::Union{Nothing, Int}=nothing
end

@inline function _parse_optional_positive_int(raw::String)::Union{Nothing, Int}
    token = strip(raw)
    isempty(token) && return nothing
    parsed = try
        parse(Int, token)
    catch
        throw(ArgumentError("Expected integer value, got '$raw'"))
    end
    if parsed <= 0
        return nothing
    end
    return parsed
end

function parse_static_vs_parallel_cli()::StaticVsParallelConfig
    profile_name = lowercase(strip(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_PROFILE", get(ENV, "SPACEAGORA_PERF_PROFILE", "full"))))
    outdir = get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_OUTDIR", STATIC_VS_PARALLEL_DEFAULT_OUTDIR)
    clean = _parse_bool_token(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_CLEAN", "1"))
    include_process = _parse_bool_token(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_INCLUDE_PROCESS", "0"))
    process_workers = _parse_optional_positive_int(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_PROCESS_WORKERS", ""))

    for arg in ARGS
        if arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(strip(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--clean=")
            clean = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--include-process=")
            include_process = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--process-workers=")
            process_workers = _parse_optional_positive_int(String(split(arg, "=", limit=2)[2]))
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..., --clean=0|1, --include-process=0|1, --process-workers=N."
            ))
        end
    end

    return StaticVsParallelConfig(
        profile=_profile_from_name(profile_name),
        outdir=abspath(outdir),
        clean=clean,
        include_process=include_process,
        process_workers=process_workers
    )
end

@inline function _default_process_workers()::Int
    available = max(1, Sys.CPU_THREADS)
    nthreads_local = max(1, Threads.nthreads())
    return max(1, fld(available, nthreads_local))
end

function _arm_specs(config::StaticVsParallelConfig)::Vector{ArmSpec}
    arms = ArmSpec[
        ArmSpec(
            mode=:serial,
            label="serial_static",
            backend="none",
            adaptive=false,
            force_inner_off=true
        ),
        ArmSpec(
            mode=:auto_static,
            label="auto_static",
            backend="auto",
            adaptive=false,
            force_inner_off=false
        ),
        ArmSpec(
            mode=:auto_adaptive,
            label="auto_adaptive",
            backend="auto",
            adaptive=true,
            force_inner_off=false
        )
    ]

    if config.include_process
        workers = isnothing(config.process_workers) ? _default_process_workers() : config.process_workers
        push!(
            arms,
            ArmSpec(
                mode=:process_static,
                label="process_static",
                backend="process",
                adaptive=false,
                force_inner_off=true,
                process_workers=workers
            )
        )
    end
    return arms
end

function _arm_env_pairs(arm::ArmSpec)::Vector{Pair{String, Union{Nothing, String}}}
    callback_mode = arm.force_inner_off ? "off" : "auto"
    pairs = Pair{String, Union{Nothing, String}}[
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => arm.backend,
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => (arm.adaptive ? "1" : "0"),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => callback_mode,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => callback_mode,
        "SPACEAGORA_MULTIBODY_PARALLEL" => callback_mode,
    ]
    if isnothing(arm.process_workers)
        push!(pairs, "SPACEAGORA_PERF_PROCS" => nothing)
    else
        push!(pairs, "SPACEAGORA_PERF_PROCS" => string(arm.process_workers))
    end
    return pairs
end

function _latest_artifact_path(outdir::String, prefix::String, profile_name::String, suffix::String)::String
    pattern = "$(prefix)_$(profile_name)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && error("No artifact found in '$outdir' for prefix='$prefix' profile='$(profile_name)' suffix='$suffix'")
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

function _arm_result_artifacts(arm::ArmSpec, config::StaticVsParallelConfig, arm_outdir::String, elapsed_s::Float64)::ModeRunArtifacts
    profile_name = config.profile.name
    raw_path = _latest_artifact_path(arm_outdir, "runtime_raw", profile_name, ".csv")
    summary_path = _latest_artifact_path(arm_outdir, "runtime_summary", profile_name, ".csv")
    orbit_raw_path = _latest_artifact_path(arm_outdir, "runtime_per_orbit_raw", profile_name, ".csv")
    orbit_summary_path = _latest_artifact_path(arm_outdir, "runtime_per_orbit_summary", profile_name, ".csv")
    report_path = _latest_artifact_path(arm_outdir, "runtime_report", profile_name, ".md")

    return ModeRunArtifacts(
        mode=arm.mode,
        backend=arm.adaptive ? string(arm.backend, "+adaptive") : string(arm.backend, "+static"),
        elapsed_s=elapsed_s,
        raw_path=raw_path,
        summary_path=summary_path,
        orbit_raw_path=orbit_raw_path,
        orbit_summary_path=orbit_summary_path,
        report_path=report_path,
        raw_df=CSV.read(raw_path, DataFrame),
        summary_df=CSV.read(summary_path, DataFrame),
        orbit_summary_df=CSV.read(orbit_summary_path, DataFrame)
    )
end

function run_arm(arm::ArmSpec, config::StaticVsParallelConfig)::ModeRunArtifacts
    arm_outdir = joinpath(config.outdir, arm.label)
    mkpath(arm_outdir)

    println()
    println("[static-vs-parallel] arm=$(arm.label) backend=$(arm.backend) adaptive=$(arm.adaptive) inner_off=$(arm.force_inner_off)")
    if !isnothing(arm.process_workers)
        println("[static-vs-parallel] arm=$(arm.label) process_workers=$(arm.process_workers)")
    end

    cmd = `$(Base.julia_cmd()) --project=$(STATIC_VS_PARALLEL_PROJECT) $(STATIC_VS_PARALLEL_RUNTIME_SCRIPT) --profile=$(config.profile.name) --outdir=$(arm_outdir)`
    env_pairs = _arm_env_pairs(arm)
    started_ns = time_ns()
    withenv(env_pairs...) do
        run(cmd)
    end
    elapsed_s = (time_ns() - started_ns) / 1e9

    artifacts = _arm_result_artifacts(arm, config, arm_outdir, elapsed_s)
    println("[static-vs-parallel] arm=$(arm.label) completed in $(round(elapsed_s; digits=3)) s")
    return artifacts
end

function _tag_arm_column(df::DataFrame, arm_label::String)::DataFrame
    tagged = copy(df)
    tagged[!, :arm] = fill(arm_label, nrow(tagged))
    return tagged
end

function _write_static_vs_parallel_report(
    path::String,
    config::StaticVsParallelConfig,
    arms::Vector{ArmSpec},
    overview_df::DataFrame,
    comparison_df::DataFrame,
    artifacts::Vector{ModeRunArtifacts};
    comparison_plot_paths::Vector{String}=String[],
    comparison_csv_path::String="",
    overview_csv_path::String="",
    merged_raw_csv_path::String="",
    merged_summary_csv_path::String=""
)
    generated = string(now(UTC))
    nthreads = Threads.nthreads()
    cpu_threads = Sys.CPU_THREADS
    arm_labels = join([arm.label for arm in arms], ", ")

    open(path, "w") do io
        println(io, "# SpaceAGORA Static vs Parallel Comparison")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Threads in wrapper process: `$nthreads`")
        println(io, "- Detected CPU threads: `$cpu_threads`")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Arms executed: `$arm_labels`")
        println(io)

        println(io, "## Arm Configuration")
        println(io)
        println(io, "| Arm | Backend | Adaptive Policy | Inner Callback/Multibody Mode | Process Workers |")
        println(io, "|---|---|---:|---|---:|")
        for arm in arms
            workers = isnothing(arm.process_workers) ? "n/a" : string(arm.process_workers)
            inner = arm.force_inner_off ? "off" : "auto"
            println(io, "| $(arm.label) | $(arm.backend) | $(arm.adaptive ? "on" : "off") | $(inner) | $(workers) |")
        end
        println(io)

        println(io, "## Mode Overview")
        println(io)
        _write_markdown_table(io, overview_df)
        println(io)

        println(io, "## Merged Comparison (Scenario Means)")
        println(io)
        println(io, "_Note: `serial_total_time_mean_s` corresponds to arm `serial_static`._")
        println(io)
        _write_markdown_table(io, comparison_df)
        println(io)

        println(io, "## Comparison Plots")
        println(io)
        println(io, "- Plot artifacts generated: `$(length(comparison_plot_paths))`")
        for plot_path in comparison_plot_paths
            println(io, "- `$(plot_path)`")
        end
        println(io)

        println(io, "## Output Files")
        println(io)
        println(io, "- Comparison CSV: `$(comparison_csv_path)`")
        println(io, "- Overview CSV: `$(overview_csv_path)`")
        println(io, "- Merged raw CSV: `$(merged_raw_csv_path)`")
        println(io, "- Merged summary CSV: `$(merged_summary_csv_path)`")
        for artifact in artifacts
            println(io, "- Arm `$(artifact.mode)` raw: `$(artifact.raw_path)`")
            println(io, "- Arm `$(artifact.mode)` summary: `$(artifact.summary_path)`")
            println(io, "- Arm `$(artifact.mode)` per-orbit summary: `$(artifact.orbit_summary_path)`")
            println(io, "- Arm `$(artifact.mode)` report: `$(artifact.report_path)`")
        end
        println(io)

        println(io, "## Reproducibility")
        println(io)
        println(io, "```bash")
        process_suffix = config.include_process ? " --include-process=1$(isnothing(config.process_workers) ? "" : " --process-workers=$(config.process_workers)")" : ""
        println(
            io,
            "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=.AGORA test/performance_static_vs_parallel.jl --profile=$(config.profile.name) --outdir=$(config.outdir) --clean=1$(process_suffix)"
        )
        println(io, "```")
    end
end

function main_static_vs_parallel()
    config = parse_static_vs_parallel_cli()
    if config.clean
        rm(STATIC_VS_PARALLEL_PERF_ROOT; recursive=true, force=true)
    end
    mkpath(config.outdir)

    arms = _arm_specs(config)
    println("Static vs parallel profile: $(config.profile.name)")
    println("Wrapper outdir: $(config.outdir)")
    println("Clean performance root: $(config.clean)")
    println("Arms: $(join([arm.label for arm in arms], ", "))")
    println("Wrapper Threads.nthreads()=$(Threads.nthreads()), Sys.CPU_THREADS=$(Sys.CPU_THREADS)")

    artifacts = ModeRunArtifacts[]
    for arm in arms
        push!(artifacts, run_arm(arm, config))
    end

    comparison_df = build_comparison_table(artifacts)
    overview_df = build_mode_overview(artifacts)
    arm_map = Dict(arm.mode => arm for arm in arms)
    overview_df.mode = [
        haskey(arm_map, Symbol(mode_name)) ? arm_map[Symbol(mode_name)].label : String(mode_name)
        for mode_name in String.(overview_df.mode)
    ]

    merged_raw_df = DataFrame()
    merged_summary_df = DataFrame()
    for artifact in artifacts
        label = arm_map[artifact.mode].label
        merged_raw_df = vcat(merged_raw_df, _tag_arm_column(artifact.raw_df, label); cols=:union)
        merged_summary_df = vcat(merged_summary_df, _tag_arm_column(artifact.summary_df, label); cols=:union)
    end

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    comparison_path = joinpath(config.outdir, "static_vs_parallel_compare_summary_$(config.profile.name)_$(stamp).csv")
    overview_path = joinpath(config.outdir, "static_vs_parallel_mode_overview_$(config.profile.name)_$(stamp).csv")
    merged_raw_path = joinpath(config.outdir, "static_vs_parallel_raw_merged_$(config.profile.name)_$(stamp).csv")
    merged_summary_path = joinpath(config.outdir, "static_vs_parallel_summary_merged_$(config.profile.name)_$(stamp).csv")
    report_path = joinpath(config.outdir, "static_vs_parallel_report_$(config.profile.name)_$(stamp).md")
    comparison_plot_paths = generate_pipeline_comparison_plots(config.outdir, config.profile, stamp, overview_df, comparison_df)

    CSV.write(comparison_path, comparison_df)
    CSV.write(overview_path, overview_df)
    CSV.write(merged_raw_path, merged_raw_df)
    CSV.write(merged_summary_path, merged_summary_df)
    _write_static_vs_parallel_report(
        report_path,
        config,
        arms,
        overview_df,
        comparison_df,
        artifacts;
        comparison_plot_paths=comparison_plot_paths,
        comparison_csv_path=comparison_path,
        overview_csv_path=overview_path,
        merged_raw_csv_path=merged_raw_path,
        merged_summary_csv_path=merged_summary_path
    )

    println()
    println("Static vs parallel comparison complete.")
    println("Comparison summary: $comparison_path")
    println("Mode overview: $overview_path")
    println("Merged raw: $merged_raw_path")
    println("Merged summary: $merged_summary_path")
    println("Comparison plots: $(length(comparison_plot_paths))")
    println("Report: $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_static_vs_parallel()
end
