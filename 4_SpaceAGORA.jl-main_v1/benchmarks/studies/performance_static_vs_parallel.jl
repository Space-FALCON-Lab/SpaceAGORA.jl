const _WRAPPER_DIR = @__DIR__
include(joinpath(_WRAPPER_DIR, "performance_paper_pipeline.jl"))

using Random

const STATIC_VS_PARALLEL_PERF_ROOT = joinpath(REPO_ROOT, "output", "performance")
const STATIC_VS_PARALLEL_DEFAULT_OUTDIR = joinpath(STATIC_VS_PARALLEL_PERF_ROOT, "static_vs_parallel")
const STATIC_VS_PARALLEL_RUNTIME_SCRIPT = joinpath(_WRAPPER_DIR, "performance_runtime_analysis.jl")
const STATIC_VS_PARALLEL_PROJECT = REPO_ROOT
const STATIC_VS_PARALLEL_POLICY_KEYS = (:outer_pinned, :full_auto)

Base.@kwdef struct StaticVsParallelConfig
    profile::ProfileSpec
    outdir::String
    clean::Bool
    include_process::Bool
    process_workers::Union{Nothing, Int}
    passes::Int
    randomize_arm_order::Bool
    random_seed::Int
    policy_matrices::Vector{Symbol}
    include_control_stress_per_orbit::Bool
    control_stress_repeats_full::Int
    control_stress_warmup_full::Int
end

Base.@kwdef struct ArmSpec
    mode::Symbol
    label::String
    backend::String
    adaptive::Bool
    process_workers::Union{Nothing, Int}=nothing
end

Base.@kwdef struct PolicyMatrixSpec
    key::Symbol
    label::String
    density_mode::String
    control_mode::String
    multibody_mode::String
    effector_mode::String
end

Base.@kwdef struct ArmPassResult
    pass::Int
    arm::ArmSpec
    artifact::ModeRunArtifacts
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

@inline function _parse_positive_int(raw::String, name::String)::Int
    token = strip(raw)
    parsed = try
        parse(Int, token)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    parsed > 0 || throw(ArgumentError("$name must be > 0, got '$parsed'"))
    return parsed
end

@inline function _parse_nonnegative_int(raw::String, name::String)::Int
    token = strip(raw)
    parsed = try
        parse(Int, token)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    parsed >= 0 || throw(ArgumentError("$name must be >= 0, got '$parsed'"))
    return parsed
end

function _parse_policy_matrices(raw::AbstractString)::Vector{Symbol}
    isempty(strip(raw)) && throw(ArgumentError("Policy matrix list cannot be empty."))
    seen = Set{Symbol}()
    parsed = Symbol[]
    for token in split(raw, ",")
        key = Symbol(lowercase(strip(token)))
        key in STATIC_VS_PARALLEL_POLICY_KEYS || throw(ArgumentError("Unsupported policy matrix '$key'. Valid: outer_pinned, full_auto."))
        if !(key in seen)
            push!(parsed, key)
            push!(seen, key)
        end
    end
    isempty(parsed) && throw(ArgumentError("No valid policy matrices were parsed from '$raw'."))
    return parsed
end

function parse_static_vs_parallel_cli()::StaticVsParallelConfig
    profile_name = lowercase(strip(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_PROFILE", get(ENV, "SPACEAGORA_PERF_PROFILE", "full"))))
    outdir = get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_OUTDIR", STATIC_VS_PARALLEL_DEFAULT_OUTDIR)
    clean = _parse_bool_token(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_CLEAN", "1"))
    include_process = _parse_bool_token(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_INCLUDE_PROCESS", "0"))
    process_workers = _parse_optional_positive_int(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_PROCESS_WORKERS", ""))
    passes = _parse_positive_int(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_PASSES", "3"), "SPACEAGORA_STATIC_VS_PARALLEL_PASSES")
    randomize_arm_order = _parse_bool_token(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_RANDOMIZE_ARM_ORDER", "1"))
    random_seed = _parse_nonnegative_int(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_SEED", "20260301"), "SPACEAGORA_STATIC_VS_PARALLEL_SEED")
    policy_matrices = _parse_policy_matrices(get(ENV, "SPACEAGORA_STATIC_VS_PARALLEL_POLICY_MATRICES", "outer_pinned,full_auto"))
    include_control_stress_per_orbit = _parse_bool_token(get(ENV, "SPACEAGORA_PERF_INCLUDE_CONTROL_STRESS_PER_ORBIT", "1"))
    control_stress_repeats_full = _parse_positive_int(get(ENV, "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL", "3"), "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL")
    control_stress_warmup_full = _parse_nonnegative_int(get(ENV, "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL", "1"), "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL")

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
        elseif startswith(arg, "--passes=")
            passes = _parse_positive_int(String(split(arg, "=", limit=2)[2]), "--passes")
        elseif startswith(arg, "--randomize-arm-order=")
            randomize_arm_order = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--seed=")
            random_seed = _parse_nonnegative_int(String(split(arg, "=", limit=2)[2]), "--seed")
        elseif startswith(arg, "--policy-matrices=")
            policy_matrices = _parse_policy_matrices(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--include-control-stress-per-orbit=")
            include_control_stress_per_orbit = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--control-stress-repeats-full=")
            control_stress_repeats_full = _parse_positive_int(String(split(arg, "=", limit=2)[2]), "--control-stress-repeats-full")
        elseif startswith(arg, "--control-stress-warmup-full=")
            control_stress_warmup_full = _parse_nonnegative_int(String(split(arg, "=", limit=2)[2]), "--control-stress-warmup-full")
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..., --clean=0|1, --include-process=0|1, --process-workers=N, --passes=N, --randomize-arm-order=0|1, --seed=N, --policy-matrices=outer_pinned,full_auto, --include-control-stress-per-orbit=0|1, --control-stress-repeats-full=N, --control-stress-warmup-full=N."
            ))
        end
    end

    return StaticVsParallelConfig(
        profile=_profile_from_name(profile_name),
        outdir=abspath(outdir),
        clean=clean,
        include_process=include_process,
        process_workers=process_workers,
        passes=passes,
        randomize_arm_order=randomize_arm_order,
        random_seed=random_seed,
        policy_matrices=policy_matrices,
        include_control_stress_per_orbit=include_control_stress_per_orbit,
        control_stress_repeats_full=control_stress_repeats_full,
        control_stress_warmup_full=control_stress_warmup_full
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
            adaptive=false
        ),
        ArmSpec(
            mode=:auto_static,
            label="auto_static",
            backend="auto",
            adaptive=false
        ),
        ArmSpec(
            mode=:auto_adaptive,
            label="auto_adaptive",
            backend="auto",
            adaptive=true
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
                process_workers=workers
            )
        )
    end
    return arms
end

@inline function _policy_matrix_spec(key::Symbol)::PolicyMatrixSpec
    if key == :outer_pinned
        return PolicyMatrixSpec(
            key=:outer_pinned,
            label="outer_pinned",
            density_mode="off",
            control_mode="off",
            multibody_mode="off",
            effector_mode="off"
        )
    elseif key == :full_auto
        return PolicyMatrixSpec(
            key=:full_auto,
            label="full_auto",
            density_mode="auto",
            control_mode="auto",
            multibody_mode="auto",
            effector_mode="auto"
        )
    end
    throw(ArgumentError("Unsupported policy matrix '$key'."))
end

function _policy_matrix_specs(config::StaticVsParallelConfig)::Vector{PolicyMatrixSpec}
    return [_policy_matrix_spec(key) for key in config.policy_matrices]
end

function _arm_env_pairs(
    arm::ArmSpec,
    matrix::PolicyMatrixSpec,
    config::StaticVsParallelConfig
)::Vector{Pair{String, Union{Nothing, String}}}
    pairs = Pair{String, Union{Nothing, String}}[
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => arm.backend,
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => (arm.adaptive ? "1" : "0"),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => matrix.density_mode,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => matrix.control_mode,
        "SPACEAGORA_MULTIBODY_PARALLEL" => matrix.multibody_mode,
        "SPACEAGORA_EFFECTOR_PARALLEL" => matrix.effector_mode,
        "SPACEAGORA_PERF_INCLUDE_CONTROL_STRESS_PER_ORBIT" => (config.include_control_stress_per_orbit ? "1" : "0"),
        "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL" => string(config.control_stress_repeats_full),
        "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL" => string(config.control_stress_warmup_full)
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

function _latest_artifact_path_optional(outdir::String, prefix::String, profile_name::String, suffix::String)::Union{Nothing, String}
    pattern = "$(prefix)_$(profile_name)_"
    candidates = String[]
    for entry in readdir(outdir)
        if startswith(entry, pattern) && endswith(entry, suffix)
            push!(candidates, joinpath(outdir, entry))
        end
    end
    isempty(candidates) && return nothing
    sort!(candidates; by=path -> stat(path).mtime)
    return last(candidates)
end

function _stage_elapsed_s(
    stage_timing_path::Union{Nothing, String},
    fallback_total_s::Float64
)::Tuple{Float64, Float64, Float64, Float64}
    if isnothing(stage_timing_path)
        return fallback_total_s, 0.0, 0.0, fallback_total_s
    end

    stage_df = CSV.read(stage_timing_path, DataFrame)
    bench_elapsed_s = 0.0
    split_gate_elapsed_s = 0.0
    orbit_elapsed_s = 0.0
    total_elapsed_s = fallback_total_s

    if :stage in names(stage_df) && :elapsed_s in names(stage_df)
        for row in eachrow(stage_df)
            stage_name = String(row.stage)
            elapsed = Float64(row.elapsed_s)
            if stage_name == "run_benchmarks"
                bench_elapsed_s = elapsed
            elseif stage_name == "run_split_rollout_gate"
                split_gate_elapsed_s = elapsed
            elseif stage_name == "run_per_orbit"
                orbit_elapsed_s = elapsed
            elseif stage_name == "total"
                total_elapsed_s = elapsed
            end
        end
    end

    if total_elapsed_s <= 0.0
        total_elapsed_s = bench_elapsed_s + split_gate_elapsed_s + orbit_elapsed_s
    end
    if total_elapsed_s <= 0.0
        total_elapsed_s = fallback_total_s
    end
    if bench_elapsed_s <= 0.0 && split_gate_elapsed_s <= 0.0 && orbit_elapsed_s <= 0.0
        bench_elapsed_s = total_elapsed_s
    end

    return bench_elapsed_s, split_gate_elapsed_s, orbit_elapsed_s, total_elapsed_s
end

function _arm_result_artifacts(
    arm::ArmSpec,
    config::StaticVsParallelConfig,
    arm_outdir::String,
    elapsed_s::Float64
)::ModeRunArtifacts
    profile_name = config.profile.name
    raw_path = _latest_artifact_path(arm_outdir, "runtime_raw", profile_name, ".csv")
    summary_path = _latest_artifact_path(arm_outdir, "runtime_summary", profile_name, ".csv")
    orbit_raw_path = _latest_artifact_path(arm_outdir, "runtime_per_orbit_raw", profile_name, ".csv")
    orbit_summary_path = _latest_artifact_path(arm_outdir, "runtime_per_orbit_summary", profile_name, ".csv")
    stage_timing_path = _latest_artifact_path_optional(arm_outdir, "runtime_stage_timing", profile_name, ".csv")
    hardware_info_path = _latest_artifact_path_optional(arm_outdir, "runtime_hardware_info", profile_name, ".csv")
    split_gate_csv_path = _latest_artifact_path_optional(arm_outdir, "split_rollout_gate", profile_name, ".csv")
    split_gate_report_path = _latest_artifact_path_optional(arm_outdir, "split_rollout_gate", profile_name, ".md")
    report_path = _latest_artifact_path(arm_outdir, "runtime_report", profile_name, ".md")

    bench_elapsed_s, split_gate_elapsed_s, orbit_elapsed_s, total_elapsed_s = _stage_elapsed_s(stage_timing_path, elapsed_s)
    split_gate_df = isnothing(split_gate_csv_path) ? nothing : CSV.read(split_gate_csv_path, DataFrame)

    return ModeRunArtifacts(
        mode=arm.mode,
        backend=arm.adaptive ? string(arm.backend, "+adaptive") : string(arm.backend, "+static"),
        elapsed_s=total_elapsed_s,
        bench_elapsed_s=bench_elapsed_s,
        orbit_elapsed_s=orbit_elapsed_s,
        raw_path=raw_path,
        summary_path=summary_path,
        orbit_raw_path=orbit_raw_path,
        orbit_summary_path=orbit_summary_path,
        report_path=report_path,
        stage_timing_path=isnothing(stage_timing_path) ? "" : stage_timing_path,
        hardware_info_path=isnothing(hardware_info_path) ? "" : hardware_info_path,
        split_gate_elapsed_s=split_gate_elapsed_s,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_report_path=split_gate_report_path,
        split_gate_df=split_gate_df,
        raw_df=CSV.read(raw_path, DataFrame),
        summary_df=CSV.read(summary_path, DataFrame),
        orbit_summary_df=CSV.read(orbit_summary_path, DataFrame)
    )
end

function run_arm(
    arm::ArmSpec,
    matrix::PolicyMatrixSpec,
    config::StaticVsParallelConfig,
    pass_idx::Int,
    pass_outdir::String
)::ModeRunArtifacts
    arm_outdir = joinpath(pass_outdir, arm.label)
    mkpath(arm_outdir)

    println()
    println("[static-vs-parallel] matrix=$(matrix.label) pass=$(pass_idx) arm=$(arm.label) backend=$(arm.backend) adaptive=$(arm.adaptive)")
    if !isnothing(arm.process_workers)
        println("[static-vs-parallel] matrix=$(matrix.label) pass=$(pass_idx) arm=$(arm.label) process_workers=$(arm.process_workers)")
    end

    cmd = `$(Base.julia_cmd()) --project=$(STATIC_VS_PARALLEL_PROJECT) $(STATIC_VS_PARALLEL_RUNTIME_SCRIPT) --profile=$(config.profile.name) --outdir=$(arm_outdir)`
    env_pairs = _arm_env_pairs(arm, matrix, config)
    started_ns = time_ns()
    withenv(env_pairs...) do
        run(cmd)
    end
    elapsed_s = (time_ns() - started_ns) / 1e9

    artifacts = _arm_result_artifacts(arm, config, arm_outdir, elapsed_s)
    println(
        "[static-vs-parallel] matrix=$(matrix.label) pass=$(pass_idx) arm=$(arm.label) completed total=$(round(artifacts.elapsed_s; digits=3)) s " *
        "(run_benchmarks=$(round(artifacts.bench_elapsed_s; digits=3)) s, split_gate=$(round(artifacts.split_gate_elapsed_s; digits=3)) s, per_orbit=$(round(artifacts.orbit_elapsed_s; digits=3)) s)"
    )
    return artifacts
end

function _tag_arm_column(
    df::DataFrame,
    arm_label::String;
    pass_idx::Union{Nothing, Int}=nothing,
    matrix_key::Union{Nothing, Symbol}=nothing
)::DataFrame
    tagged = copy(df)
    tagged[!, :arm] = fill(arm_label, nrow(tagged))
    if !isnothing(pass_idx)
        tagged[!, :pass] = fill(pass_idx, nrow(tagged))
    end
    if !isnothing(matrix_key)
        tagged[!, :policy_matrix] = fill(String(matrix_key), nrow(tagged))
    end
    return tagged
end

function _ordered_arms(
    arms::Vector{ArmSpec},
    rng::Random.AbstractRNG,
    randomize::Bool
)::Vector{ArmSpec}
    if !randomize || length(arms) <= 1
        return copy(arms)
    end
    perm = randperm(rng, length(arms))
    return arms[perm]
end

function _write_aggregate_arm_report(
    path::String,
    config::StaticVsParallelConfig,
    matrix::PolicyMatrixSpec,
    arm::ArmSpec,
    runs::Vector{ArmPassResult},
    bench_elapsed_s::Float64,
    orbit_elapsed_s::Float64,
    total_elapsed_s::Float64;
    split_gate_csv_path::Union{Nothing, String}=nothing,
    split_gate_pass_rows::Int=0,
    split_gate_total_rows::Int=0
)
    open(path, "w") do io
        println(io, "# Static vs Parallel Aggregated Arm Report")
        println(io)
        println(io, "- Generated (UTC): $(now(UTC))")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Policy matrix: `$(matrix.key)`")
        println(io, "- Arm: `$(arm.label)`")
        println(io, "- Aggregated passes: `$(length(runs))`")
        println(io, "- Mean run_benchmarks elapsed: `$(round(bench_elapsed_s; digits=3)) s`")
        println(io, "- Mean per-orbit elapsed: `$(round(orbit_elapsed_s; digits=3)) s`")
        println(io, "- Mean total elapsed: `$(round(total_elapsed_s; digits=3)) s`")
        if split_gate_total_rows > 0
            println(io, "- Split rollout gate pass rows: `$(split_gate_pass_rows)/$(split_gate_total_rows)`")
        else
            println(io, "- Split rollout gate rows: `0`")
        end
        if !(split_gate_csv_path === nothing)
            println(io, "- Split rollout gate aggregated CSV: `$(split_gate_csv_path)`")
        end
        println(io)
        println(io, "## Source Runs")
        println(io)
        for run in runs
            println(io, "- pass=$(run.pass): raw=`$(run.artifact.raw_path)`, per-orbit raw=`$(run.artifact.orbit_raw_path)`, report=`$(run.artifact.report_path)`")
        end
    end
end

function _aggregate_arm_artifacts(
    matrix_outdir::String,
    matrix::PolicyMatrixSpec,
    config::StaticVsParallelConfig,
    arm::ArmSpec,
    runs::Vector{ArmPassResult}
)::ModeRunArtifacts
    isempty(runs) && error("No pass data collected for arm '$(arm.label)'.")

    raw_df = DataFrame()
    orbit_raw_df = DataFrame()
    split_gate_df = DataFrame()
    for run in runs
        raw_df = vcat(raw_df, _tag_arm_column(run.artifact.raw_df, arm.label; pass_idx=run.pass, matrix_key=matrix.key); cols=:union)
        local_orbit_raw = CSV.read(run.artifact.orbit_raw_path, DataFrame)
        orbit_raw_df = vcat(orbit_raw_df, _tag_arm_column(local_orbit_raw, arm.label; pass_idx=run.pass, matrix_key=matrix.key); cols=:union)
        if !(run.artifact.split_gate_df === nothing)
            split_gate_df = vcat(
                split_gate_df,
                _tag_arm_column(run.artifact.split_gate_df, arm.label; pass_idx=run.pass, matrix_key=matrix.key);
                cols=:union
            )
        end
    end

    summary_df = summarize_results(raw_df)
    orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)

    bench_elapsed_s = mean([run.artifact.bench_elapsed_s for run in runs])
    split_gate_elapsed_s = mean([run.artifact.split_gate_elapsed_s for run in runs])
    orbit_elapsed_s = mean([run.artifact.orbit_elapsed_s for run in runs])
    total_elapsed_s = bench_elapsed_s + split_gate_elapsed_s + orbit_elapsed_s

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    agg_outdir = joinpath(matrix_outdir, "aggregate", arm.label)
    mkpath(agg_outdir)

    raw_path = joinpath(agg_outdir, "runtime_raw_agg_$(config.profile.name)_$(arm.label)_$(stamp).csv")
    summary_path = joinpath(agg_outdir, "runtime_summary_agg_$(config.profile.name)_$(arm.label)_$(stamp).csv")
    orbit_raw_path = joinpath(agg_outdir, "runtime_per_orbit_raw_agg_$(config.profile.name)_$(arm.label)_$(stamp).csv")
    orbit_summary_path = joinpath(agg_outdir, "runtime_per_orbit_summary_agg_$(config.profile.name)_$(arm.label)_$(stamp).csv")
    stage_timing_path = joinpath(agg_outdir, "runtime_stage_timing_agg_$(config.profile.name)_$(arm.label)_$(stamp).csv")
    hardware_info_path = joinpath(agg_outdir, "runtime_hardware_info_agg_$(config.profile.name)_$(arm.label)_$(stamp).csv")
    split_gate_csv_path = nrow(split_gate_df) > 0 ? joinpath(agg_outdir, "split_rollout_gate_agg_$(config.profile.name)_$(arm.label)_$(stamp).csv") : nothing
    report_path = joinpath(agg_outdir, "runtime_report_agg_$(config.profile.name)_$(arm.label)_$(stamp).md")
    split_gate_pass_rows = (nrow(split_gate_df) > 0 && (:pass_all in names(split_gate_df))) ? count(Bool.(split_gate_df.pass_all)) : 0
    split_gate_total_rows = nrow(split_gate_df)
    hw = _runtime_hardware_snapshot()
    hardware_info_df = DataFrame([
        (
            profile=config.profile.name,
            policy_matrix=String(matrix.key),
            arm=arm.label,
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
    stage_timing_df = DataFrame(
        stage=["run_benchmarks", "run_split_rollout_gate", "run_per_orbit", "total"],
        elapsed_s=[bench_elapsed_s, split_gate_elapsed_s, orbit_elapsed_s, total_elapsed_s]
    )

    CSV.write(raw_path, raw_df)
    CSV.write(summary_path, summary_df)
    CSV.write(orbit_raw_path, orbit_raw_df)
    CSV.write(orbit_summary_path, orbit_summary_df)
    CSV.write(stage_timing_path, stage_timing_df)
    CSV.write(hardware_info_path, hardware_info_df)
    if !(split_gate_csv_path === nothing)
        CSV.write(split_gate_csv_path, split_gate_df)
    end
    _write_aggregate_arm_report(
        report_path,
        config,
        matrix,
        arm,
        runs,
        bench_elapsed_s,
        orbit_elapsed_s,
        total_elapsed_s;
        split_gate_csv_path=split_gate_csv_path,
        split_gate_pass_rows=split_gate_pass_rows,
        split_gate_total_rows=split_gate_total_rows
    )

    backend = runs[1].artifact.backend
    return ModeRunArtifacts(
        mode=arm.mode,
        backend=backend,
        elapsed_s=total_elapsed_s,
        bench_elapsed_s=bench_elapsed_s,
        orbit_elapsed_s=orbit_elapsed_s,
        raw_path=raw_path,
        summary_path=summary_path,
        orbit_raw_path=orbit_raw_path,
        orbit_summary_path=orbit_summary_path,
        report_path=report_path,
        stage_timing_path=stage_timing_path,
        hardware_info_path=hardware_info_path,
        split_gate_elapsed_s=split_gate_elapsed_s,
        split_gate_csv_path=split_gate_csv_path,
        split_gate_df=(split_gate_total_rows > 0 ? split_gate_df : nothing),
        raw_df=raw_df,
        summary_df=summary_df,
        orbit_summary_df=orbit_summary_df
    )
end

function _write_static_vs_parallel_report(
    path::String,
    config::StaticVsParallelConfig,
    matrix::PolicyMatrixSpec,
    arms::Vector{ArmSpec},
    overview_df::DataFrame,
    comparison_df::DataFrame,
    artifacts::Vector{ModeRunArtifacts},
    run_order_df::DataFrame;
    comparison_plot_paths::Vector{String}=String[],
    comparison_csv_path::String="",
    overview_csv_path::String="",
    merged_raw_csv_path::String="",
    merged_summary_csv_path::String="",
    run_order_csv_path::String=""
)
    generated = string(now(UTC))
    nthreads = Threads.nthreads()
    cpu_threads = Sys.CPU_THREADS
    arm_labels = join([arm.label for arm in arms], ", ")
    matrix_scenarios = Set(String.(comparison_df.scenario))
    raw_scenarios = Set{String}()
    for artifact in artifacts
        for s in String.(artifact.summary_df.scenario)
            push!(raw_scenarios, s)
        end
    end
    missing_from_matrix = sort(collect(setdiff(raw_scenarios, matrix_scenarios)))
    split_gate_total_rows = 0
    split_gate_pass_rows = 0
    for artifact in artifacts
        if !(artifact.split_gate_df === nothing) && (:pass_all in names(artifact.split_gate_df))
            split_gate_total_rows += nrow(artifact.split_gate_df)
            split_gate_pass_rows += count(Bool.(artifact.split_gate_df.pass_all))
        end
    end

    open(path, "w") do io
        println(io, "# SpaceAGORA Static vs Parallel Comparison")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Threads in wrapper process: `$nthreads`")
        println(io, "- Detected CPU threads: `$cpu_threads`")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Policy matrix: `$(matrix.key)`")
        println(io, "- Arms executed: `$arm_labels`")
        println(io, "- Passes: `$(config.passes)`")
        println(io, "- Randomized arm order: `$(config.randomize_arm_order)` (seed=`$(config.random_seed)`)")
        println(io, "- Include control_stress in per-orbit stage: `$(config.include_control_stress_per_orbit)`")
        println(io, "- Control_stress full schedule override: warmup=`$(config.control_stress_warmup_full)`, repeats=`$(config.control_stress_repeats_full)`")
        println(io)

        println(io, "## Policy Configuration")
        println(io)
        println(io, "| Policy Matrix | Density Callback | Control Callback | Multibody | Dynamic Effector |")
        println(io, "|---|---|---|---|---|")
        println(io, "| $(matrix.key) | $(matrix.density_mode) | $(matrix.control_mode) | $(matrix.multibody_mode) | $(matrix.effector_mode) |")
        println(io)

        println(io, "## Arm Configuration")
        println(io)
        println(io, "| Arm | Backend | Adaptive Policy | Process Workers |")
        println(io, "|---|---|---:|---:|")
        for arm in arms
            workers = isnothing(arm.process_workers) ? "n/a" : string(arm.process_workers)
            println(io, "| $(arm.label) | $(arm.backend) | $(arm.adaptive ? "on" : "off") | $(workers) |")
        end
        println(io)

        println(io, "## Hardware Snapshot")
        println(io)
        println(io, "| Arm | Machine Label | Hardware Class | CPU Threads | Julia Threads |")
        println(io, "|---|---|---|---:|---:|")
        for artifact in artifacts
            raw = artifact.raw_df
            machine = (:machine_label in names(raw) && nrow(raw) > 0) ? string(raw.machine_label[1]) : "n/a"
            hw_class = (:hardware_class in names(raw) && nrow(raw) > 0) ? string(raw.hardware_class[1]) : "n/a"
            cpu_t = (:cpu_threads in names(raw) && nrow(raw) > 0) ? string(raw.cpu_threads[1]) : "n/a"
            julia_t = (:julia_threads in names(raw) && nrow(raw) > 0) ? string(raw.julia_threads[1]) : "n/a"
            println(io, "| $(artifact.mode) | $(machine) | $(hw_class) | $(cpu_t) | $(julia_t) |")
        end
        println(io)

        println(io, "## Scenario Coverage")
        println(io)
        if isempty(missing_from_matrix)
            println(io, "- No scenarios are missing from the comparison matrix.")
        else
            println(io, "- Scenarios missing from comparison matrix: `$(join(missing_from_matrix, ", "))`")
        end
        println(io)

        println(io, "## Run Order By Pass")
        println(io)
        _write_markdown_table(io, run_order_df)
        println(io)

        println(io, "## Mode Overview (Aggregated Across Passes)")
        println(io)
        _write_markdown_table(io, overview_df)
        println(io)

        println(io, "## Merged Comparison (Scenario Means + CI)")
        println(io)
        println(io, "_Note: `serial_total_time_mean_s` corresponds to arm `serial_static`._")
        println(io)
        _write_markdown_table(io, comparison_df)
        println(io)

        println(io, "## Fidelity Guardrails")
        println(io)
        if split_gate_total_rows > 0
            pass_rate = 100.0 * split_gate_pass_rows / split_gate_total_rows
            println(io, "- Split rollout gate rows (aggregated): `$(split_gate_pass_rows)/$(split_gate_total_rows)` pass (`$(round(pass_rate; digits=2))%`).")
        else
            println(io, "- Split rollout gate rows: none (gate disabled or no gate artifacts found).")
        end
        for artifact in artifacts
            if !(artifact.split_gate_csv_path === nothing)
                println(io, "- Arm `$(artifact.mode)` split gate CSV: `$(artifact.split_gate_csv_path)`")
            end
            if !(artifact.split_gate_report_path === nothing)
                println(io, "- Arm `$(artifact.mode)` split gate report: `$(artifact.split_gate_report_path)`")
            end
        end
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
        println(io, "- Merged raw CSV (all passes): `$(merged_raw_csv_path)`")
        println(io, "- Merged summary CSV (all passes): `$(merged_summary_csv_path)`")
        println(io, "- Run order CSV: `$(run_order_csv_path)`")
        for artifact in artifacts
            println(io, "- Arm `$(artifact.mode)` aggregated raw: `$(artifact.raw_path)`")
            println(io, "- Arm `$(artifact.mode)` aggregated summary: `$(artifact.summary_path)`")
            println(io, "- Arm `$(artifact.mode)` aggregated per-orbit summary: `$(artifact.orbit_summary_path)`")
            if !isempty(artifact.stage_timing_path)
                println(io, "- Arm `$(artifact.mode)` aggregated stage timing: `$(artifact.stage_timing_path)`")
            end
            if !isempty(artifact.hardware_info_path)
                println(io, "- Arm `$(artifact.mode)` aggregated hardware info: `$(artifact.hardware_info_path)`")
            end
            if !(artifact.split_gate_csv_path === nothing)
                println(io, "- Arm `$(artifact.mode)` aggregated split gate CSV: `$(artifact.split_gate_csv_path)`")
            end
            println(io, "- Arm `$(artifact.mode)` aggregated report: `$(artifact.report_path)`")
        end
        println(io)

        println(io, "## Reproducibility")
        println(io)
        println(io, "```bash")
        process_suffix = config.include_process ? " --include-process=1$(isnothing(config.process_workers) ? "" : " --process-workers=$(config.process_workers)")" : ""
        println(
            io,
            "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=. benchmarks/studies/performance_static_vs_parallel.jl --profile=$(config.profile.name) --outdir=$(config.outdir) --clean=1 --passes=$(config.passes) --randomize-arm-order=$(config.randomize_arm_order ? 1 : 0) --seed=$(config.random_seed) --policy-matrices=$(join(string.(config.policy_matrices), ',')) --include-control-stress-per-orbit=$(config.include_control_stress_per_orbit ? 1 : 0) --control-stress-repeats-full=$(config.control_stress_repeats_full) --control-stress-warmup-full=$(config.control_stress_warmup_full)$(process_suffix)"
        )
        println(io, "```")
    end
end

function _run_policy_matrix(
    matrix::PolicyMatrixSpec,
    config::StaticVsParallelConfig,
    arms::Vector{ArmSpec}
)
    matrix_outdir = joinpath(config.outdir, String(matrix.key))
    mkpath(matrix_outdir)

    matrix_seed = config.random_seed + Int(sum(codeunits(String(matrix.key))))
    rng = MersenneTwister(matrix_seed)

    println()
    println("[static-vs-parallel] policy matrix=$(matrix.key) outdir=$(matrix_outdir)")
    println(
        "[static-vs-parallel] policy matrix=$(matrix.key) callbacks=(density=$(matrix.density_mode), control=$(matrix.control_mode), multibody=$(matrix.multibody_mode)), " *
        "effector=$(matrix.effector_mode), passes=$(config.passes), randomize=$(config.randomize_arm_order), seed=$(matrix_seed)"
    )

    pass_results = ArmPassResult[]
    order_rows = NamedTuple[]

    for pass_idx in 1:config.passes
        pass_outdir = joinpath(matrix_outdir, "pass_$(lpad(pass_idx, 2, '0'))")
        mkpath(pass_outdir)

        ordered_arms = _ordered_arms(arms, rng, config.randomize_arm_order)
        order_line = join([arm.label for arm in ordered_arms], " -> ")
        println("[static-vs-parallel] policy matrix=$(matrix.key) pass=$(pass_idx) arm-order=$(order_line)")

        for (order_idx, arm) in enumerate(ordered_arms)
            push!(order_rows, (
                pass=pass_idx,
                order_index=order_idx,
                arm=arm.label,
                backend=arm.backend,
                adaptive=arm.adaptive
            ))
            artifact = run_arm(arm, matrix, config, pass_idx, pass_outdir)
            push!(pass_results, ArmPassResult(pass=pass_idx, arm=arm, artifact=artifact))
        end
    end

    aggregated_artifacts = ModeRunArtifacts[]
    for arm in arms
        runs = [result for result in pass_results if result.arm.label == arm.label]
        push!(aggregated_artifacts, _aggregate_arm_artifacts(matrix_outdir, matrix, config, arm, runs))
    end

    comparison_df = build_comparison_table(aggregated_artifacts)
    overview_df = build_mode_overview(aggregated_artifacts)
    arm_map = Dict(arm.mode => arm for arm in arms)
    overview_df.mode = [
        haskey(arm_map, Symbol(mode_name)) ? arm_map[Symbol(mode_name)].label : String(mode_name)
        for mode_name in String.(overview_df.mode)
    ]

    merged_raw_df = DataFrame()
    merged_summary_df = DataFrame()
    for result in pass_results
        merged_raw_df = vcat(
            merged_raw_df,
            _tag_arm_column(result.artifact.raw_df, result.arm.label; pass_idx=result.pass, matrix_key=matrix.key);
            cols=:union
        )
        merged_summary_df = vcat(
            merged_summary_df,
            _tag_arm_column(result.artifact.summary_df, result.arm.label; pass_idx=result.pass, matrix_key=matrix.key);
            cols=:union
        )
    end

    run_order_df = DataFrame(order_rows)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    comparison_path = joinpath(matrix_outdir, "static_vs_parallel_compare_summary_$(config.profile.name)_$(matrix.key)_$(stamp).csv")
    overview_path = joinpath(matrix_outdir, "static_vs_parallel_mode_overview_$(config.profile.name)_$(matrix.key)_$(stamp).csv")
    merged_raw_path = joinpath(matrix_outdir, "static_vs_parallel_raw_merged_$(config.profile.name)_$(matrix.key)_$(stamp).csv")
    merged_summary_path = joinpath(matrix_outdir, "static_vs_parallel_summary_merged_$(config.profile.name)_$(matrix.key)_$(stamp).csv")
    run_order_path = joinpath(matrix_outdir, "static_vs_parallel_run_order_$(config.profile.name)_$(matrix.key)_$(stamp).csv")
    report_path = joinpath(matrix_outdir, "static_vs_parallel_report_$(config.profile.name)_$(matrix.key)_$(stamp).md")
    comparison_plot_paths = generate_pipeline_comparison_plots(matrix_outdir, config.profile, "$(matrix.key)_$(stamp)", overview_df, comparison_df)

    CSV.write(comparison_path, comparison_df)
    CSV.write(overview_path, overview_df)
    CSV.write(merged_raw_path, merged_raw_df)
    CSV.write(merged_summary_path, merged_summary_df)
    CSV.write(run_order_path, run_order_df)
    _write_static_vs_parallel_report(
        report_path,
        config,
        matrix,
        arms,
        overview_df,
        comparison_df,
        aggregated_artifacts,
        run_order_df;
        comparison_plot_paths=comparison_plot_paths,
        comparison_csv_path=comparison_path,
        overview_csv_path=overview_path,
        merged_raw_csv_path=merged_raw_path,
        merged_summary_csv_path=merged_summary_path,
        run_order_csv_path=run_order_path
    )

    println("[static-vs-parallel] policy matrix=$(matrix.key) comparison summary: $comparison_path")
    println("[static-vs-parallel] policy matrix=$(matrix.key) mode overview: $overview_path")
    println("[static-vs-parallel] policy matrix=$(matrix.key) merged raw: $merged_raw_path")
    println("[static-vs-parallel] policy matrix=$(matrix.key) merged summary: $merged_summary_path")
    println("[static-vs-parallel] policy matrix=$(matrix.key) run order: $run_order_path")
    println("[static-vs-parallel] policy matrix=$(matrix.key) comparison plots: $(length(comparison_plot_paths))")
    println("[static-vs-parallel] policy matrix=$(matrix.key) report: $report_path")

    return (
        key=matrix.key,
        outdir=matrix_outdir,
        comparison_path=comparison_path,
        overview_path=overview_path,
        merged_raw_path=merged_raw_path,
        merged_summary_path=merged_summary_path,
        run_order_path=run_order_path,
        report_path=report_path,
        plot_count=length(comparison_plot_paths)
    )
end

function main_static_vs_parallel()
    config = parse_static_vs_parallel_cli()
    if config.clean
        rm(STATIC_VS_PARALLEL_PERF_ROOT; recursive=true, force=true)
    end
    mkpath(config.outdir)

    arms = _arm_specs(config)
    matrices = _policy_matrix_specs(config)
    println("Static vs parallel profile: $(config.profile.name)")
    println("Wrapper outdir: $(config.outdir)")
    println("Clean performance root: $(config.clean)")
    println("Arms: $(join([arm.label for arm in arms], ", "))")
    println("Policy matrices: $(join(string.(config.policy_matrices), ", "))")
    println("Passes: $(config.passes)")
    println("Randomized arm order: $(config.randomize_arm_order) (seed=$(config.random_seed))")
    println("Include control_stress in per-orbit stage: $(config.include_control_stress_per_orbit)")
    println("Control_stress full schedule override: warmup=$(config.control_stress_warmup_full), repeats=$(config.control_stress_repeats_full)")
    println("Wrapper Threads.nthreads()=$(Threads.nthreads()), Sys.CPU_THREADS=$(Sys.CPU_THREADS)")

    matrix_results = NamedTuple[]
    for matrix in matrices
        push!(matrix_results, _run_policy_matrix(matrix, config, arms))
    end

    println()
    println("Static vs parallel comparison complete.")
    for result in matrix_results
        println("Policy matrix $(result.key):")
        println("  comparison summary: $(result.comparison_path)")
        println("  mode overview: $(result.overview_path)")
        println("  merged raw: $(result.merged_raw_path)")
        println("  merged summary: $(result.merged_summary_path)")
        println("  run order: $(result.run_order_path)")
        println("  comparison plots: $(result.plot_count)")
        println("  report: $(result.report_path)")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_static_vs_parallel()
end
