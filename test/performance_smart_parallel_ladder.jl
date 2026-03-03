const _SMART_LADDER_DIR = @__DIR__
include(joinpath(_SMART_LADDER_DIR, "performance_paper_pipeline.jl"))

using Random

const SMART_LADDER_DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "performance", "smart_parallel_ladder")
const SMART_LADDER_RUNTIME_SCRIPT = joinpath(_SMART_LADDER_DIR, "performance_runtime_analysis.jl")
const SMART_LADDER_PROJECT = joinpath(REPO_ROOT, ".AGORA")

Base.@kwdef struct SmartLadderConfig
    profile::ProfileSpec
    outdir::String
    clean::Bool
    passes::Int
    randomize_rung_order::Bool
    random_seed::Int
    outer_only_backend::String
    process_workers::Union{Nothing, Int}
    include_control_stress_per_orbit::Bool
    control_stress_repeats_full::Int
    control_stress_warmup_full::Int
end

Base.@kwdef struct LadderRungSpec
    mode::Symbol
    label::String
    description::String
    matrix::Symbol
    backend::String
    inner_adaptive::Bool
    outer_route_adaptive::Bool
end

Base.@kwdef struct LadderPassResult
    pass::Int
    rung::LadderRungSpec
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
    parsed > 0 || return nothing
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

@inline function _parse_outer_only_backend(raw::AbstractString)::String
    backend = lowercase(strip(String(raw)))
    backend in ("threads", "process", "auto") || throw(
        ArgumentError("outer-only backend must be one of: threads, process, auto (got '$raw').")
    )
    return backend
end

function parse_smart_ladder_cli()::SmartLadderConfig
    profile_name = lowercase(strip(get(ENV, "SPACEAGORA_SMART_LADDER_PROFILE", get(ENV, "SPACEAGORA_PERF_PROFILE", "full"))))
    outdir = get(ENV, "SPACEAGORA_SMART_LADDER_OUTDIR", SMART_LADDER_DEFAULT_OUTDIR)
    clean = _parse_bool_token(get(ENV, "SPACEAGORA_SMART_LADDER_CLEAN", "1"))
    passes = _parse_positive_int(get(ENV, "SPACEAGORA_SMART_LADDER_PASSES", "3"), "SPACEAGORA_SMART_LADDER_PASSES")
    randomize_rung_order = _parse_bool_token(get(ENV, "SPACEAGORA_SMART_LADDER_RANDOMIZE_ORDER", "1"))
    random_seed = _parse_nonnegative_int(get(ENV, "SPACEAGORA_SMART_LADDER_SEED", "20260302"), "SPACEAGORA_SMART_LADDER_SEED")
    outer_only_backend = _parse_outer_only_backend(get(ENV, "SPACEAGORA_SMART_LADDER_OUTER_ONLY_BACKEND", "threads"))
    process_workers = _parse_optional_positive_int(get(ENV, "SPACEAGORA_SMART_LADDER_PROCESS_WORKERS", ""))
    include_control_stress_per_orbit = _parse_bool_token(get(ENV, "SPACEAGORA_PERF_INCLUDE_CONTROL_STRESS_PER_ORBIT", "1"))
    control_stress_repeats_full = _parse_positive_int(
        get(ENV, "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL", "3"),
        "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL"
    )
    control_stress_warmup_full = _parse_nonnegative_int(
        get(ENV, "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL", "1"),
        "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL"
    )

    for arg in ARGS
        if arg in ("quick", "full")
            profile_name = arg
        elseif startswith(arg, "--profile=")
            profile_name = lowercase(strip(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--clean=")
            clean = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--passes=")
            passes = _parse_positive_int(String(split(arg, "=", limit=2)[2]), "--passes")
        elseif startswith(arg, "--randomize-rung-order=")
            randomize_rung_order = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--randomize-order=")
            randomize_rung_order = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--seed=")
            random_seed = _parse_nonnegative_int(String(split(arg, "=", limit=2)[2]), "--seed")
        elseif startswith(arg, "--outer-only-backend=")
            outer_only_backend = _parse_outer_only_backend(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--process-workers=")
            process_workers = _parse_optional_positive_int(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--include-control-stress-per-orbit=")
            include_control_stress_per_orbit = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--control-stress-repeats-full=")
            control_stress_repeats_full = _parse_positive_int(
                String(split(arg, "=", limit=2)[2]),
                "--control-stress-repeats-full"
            )
        elseif startswith(arg, "--control-stress-warmup-full=")
            control_stress_warmup_full = _parse_nonnegative_int(
                String(split(arg, "=", limit=2)[2]),
                "--control-stress-warmup-full"
            )
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: [quick|full], --profile=..., --outdir=..., --clean=0|1, " *
                "--passes=N, --randomize-rung-order=0|1, --seed=N, --outer-only-backend=threads|process|auto, " *
                "--process-workers=N, --include-control-stress-per-orbit=0|1, --control-stress-repeats-full=N, " *
                "--control-stress-warmup-full=N."
            ))
        end
    end

    return SmartLadderConfig(
        profile=_profile_from_name(profile_name),
        outdir=abspath(outdir),
        clean=clean,
        passes=passes,
        randomize_rung_order=randomize_rung_order,
        random_seed=random_seed,
        outer_only_backend=outer_only_backend,
        process_workers=process_workers,
        include_control_stress_per_orbit=include_control_stress_per_orbit,
        control_stress_repeats_full=control_stress_repeats_full,
        control_stress_warmup_full=control_stress_warmup_full
    )
end

@inline function _ladder_matrix_modes(matrix::Symbol)::NamedTuple
    if matrix == :outer_pinned
        return (density="off", control="off", multibody="off", effector="off")
    elseif matrix == :full_auto
        return (density="auto", control="auto", multibody="auto", effector="auto")
    end
    throw(ArgumentError("Unsupported ladder matrix '$matrix'. Use: outer_pinned, full_auto."))
end

function _ladder_rungs(config::SmartLadderConfig)::Vector{LadderRungSpec}
    return LadderRungSpec[
        LadderRungSpec(
            mode=:serial,
            label="r0_true_serial",
            description="No outer parallelism, no inner callback parallelism, no adaptive policies.",
            matrix=:outer_pinned,
            backend="none",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_only_process,
            label="r1_b_outer_only_process",
            description="Outer-loop process parallel only; inner callback parallelism pinned off.",
            matrix=:outer_pinned,
            backend="process",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_only,
            label="r1_a_outer_only",
            description="Outer-loop parallel only; inner callback parallelism pinned off.",
            matrix=:outer_pinned,
            backend=config.outer_only_backend,
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:inner_only,
            label="r2_inner_only",
            description="Inner callback auto policies only; outer loop serial.",
            matrix=:full_auto,
            backend="none",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_inner_static,
            label="r3_outer_inner_static",
            description="Outer + inner enabled with static policy behavior (adaptive off).",
            matrix=:full_auto,
            backend="auto",
            inner_adaptive=false,
            outer_route_adaptive=false
        ),
        LadderRungSpec(
            mode=:outer_inner_adaptive,
            label="r4_outer_inner_adaptive",
            description="Outer + inner enabled with adaptive policy behavior.",
            matrix=:full_auto,
            backend="auto",
            inner_adaptive=true,
            outer_route_adaptive=true
        )
    ]
end

function _ladder_env_pairs(
    rung::LadderRungSpec,
    config::SmartLadderConfig
)::Vector{Pair{String, Union{Nothing, String}}}
    matrix = _ladder_matrix_modes(rung.matrix)
    pairs = Pair{String, Union{Nothing, String}}[
        "SPACEAGORA_PERF_PARALLEL_BACKEND" => rung.backend,
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => (rung.inner_adaptive ? "1" : "0"),
        "SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE" => (rung.outer_route_adaptive ? "1" : "0"),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => matrix.density,
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => matrix.control,
        "SPACEAGORA_MULTIBODY_PARALLEL" => matrix.multibody,
        "SPACEAGORA_EFFECTOR_PARALLEL" => matrix.effector,
        "SPACEAGORA_PERF_INCLUDE_CONTROL_STRESS_PER_ORBIT" => (config.include_control_stress_per_orbit ? "1" : "0"),
        "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL" => string(config.control_stress_repeats_full),
        "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL" => string(config.control_stress_warmup_full)
    ]
    if isnothing(config.process_workers)
        push!(pairs, "SPACEAGORA_PERF_PROCS" => nothing)
    else
        push!(pairs, "SPACEAGORA_PERF_PROCS" => string(config.process_workers))
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

@inline function _has_column(df::DataFrame, column::Symbol)::Bool
    return column in propertynames(df)
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

    if _has_column(stage_df, :stage) && _has_column(stage_df, :elapsed_s)
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

function _ladder_artifacts(
    rung::LadderRungSpec,
    config::SmartLadderConfig,
    rung_outdir::String,
    elapsed_s::Float64
)::ModeRunArtifacts
    profile_name = config.profile.name
    raw_path = _latest_artifact_path(rung_outdir, "runtime_raw", profile_name, ".csv")
    summary_path = _latest_artifact_path(rung_outdir, "runtime_summary", profile_name, ".csv")
    orbit_raw_path = _latest_artifact_path(rung_outdir, "runtime_per_orbit_raw", profile_name, ".csv")
    orbit_summary_path = _latest_artifact_path(rung_outdir, "runtime_per_orbit_summary", profile_name, ".csv")
    stage_timing_path = _latest_artifact_path_optional(rung_outdir, "runtime_stage_timing", profile_name, ".csv")
    hardware_info_path = _latest_artifact_path_optional(rung_outdir, "runtime_hardware_info", profile_name, ".csv")
    split_gate_csv_path = _latest_artifact_path_optional(rung_outdir, "split_rollout_gate", profile_name, ".csv")
    split_gate_report_path = _latest_artifact_path_optional(rung_outdir, "split_rollout_gate", profile_name, ".md")
    report_path = _latest_artifact_path(rung_outdir, "runtime_report", profile_name, ".md")

    bench_elapsed_s, split_gate_elapsed_s, orbit_elapsed_s, total_elapsed_s = _stage_elapsed_s(stage_timing_path, elapsed_s)
    split_gate_df = isnothing(split_gate_csv_path) ? nothing : CSV.read(split_gate_csv_path, DataFrame)
    backend_label =
        "$(rung.backend)|matrix=$(rung.matrix)|inner=$(rung.inner_adaptive ? "adaptive" : "static")|outer=$(rung.outer_route_adaptive ? "adaptive" : "static")"

    return ModeRunArtifacts(
        mode=rung.mode,
        backend=backend_label,
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

function run_rung(
    rung::LadderRungSpec,
    config::SmartLadderConfig,
    pass_idx::Int,
    pass_outdir::String
)::ModeRunArtifacts
    rung_outdir = joinpath(pass_outdir, rung.label)
    mkpath(rung_outdir)

    println()
    println(
        "[smart-ladder] pass=$(pass_idx) rung=$(rung.label) backend=$(rung.backend) matrix=$(rung.matrix) " *
        "inner_adaptive=$(rung.inner_adaptive) outer_route_adaptive=$(rung.outer_route_adaptive)"
    )
    cmd = `$(Base.julia_cmd()) --project=$(SMART_LADDER_PROJECT) $(SMART_LADDER_RUNTIME_SCRIPT) --profile=$(config.profile.name) --outdir=$(rung_outdir)`
    env_pairs = _ladder_env_pairs(rung, config)

    started_ns = time_ns()
    withenv(env_pairs...) do
        run(cmd)
    end
    elapsed_s = (time_ns() - started_ns) / 1e9

    artifacts = _ladder_artifacts(rung, config, rung_outdir, elapsed_s)
    println(
        "[smart-ladder] pass=$(pass_idx) rung=$(rung.label) completed total=$(round(artifacts.elapsed_s; digits=3)) s " *
        "(run_benchmarks=$(round(artifacts.bench_elapsed_s; digits=3)) s, split_gate=$(round(artifacts.split_gate_elapsed_s; digits=3)) s, " *
        "per_orbit=$(round(artifacts.orbit_elapsed_s; digits=3)) s)"
    )
    return artifacts
end

function _tag_rung_column(
    df::DataFrame,
    rung::LadderRungSpec;
    pass_idx::Union{Nothing, Int}=nothing
)::DataFrame
    tagged = copy(df)
    tagged[!, :rung] = fill(rung.label, nrow(tagged))
    tagged[!, :rung_mode] = fill(String(rung.mode), nrow(tagged))
    tagged[!, :rung_matrix] = fill(String(rung.matrix), nrow(tagged))
    tagged[!, :rung_backend] = fill(rung.backend, nrow(tagged))
    tagged[!, :rung_inner_adaptive] = fill(rung.inner_adaptive, nrow(tagged))
    tagged[!, :rung_outer_route_adaptive] = fill(rung.outer_route_adaptive, nrow(tagged))
    if !isnothing(pass_idx)
        tagged[!, :pass] = fill(pass_idx, nrow(tagged))
    end
    return tagged
end

function _ordered_rungs(
    rungs::Vector{LadderRungSpec},
    rng::Random.AbstractRNG,
    randomize::Bool
)::Vector{LadderRungSpec}
    if !randomize || length(rungs) <= 1
        return copy(rungs)
    end
    perm = randperm(rng, length(rungs))
    return rungs[perm]
end

@inline function _safe_ratio(num::Real, den::Real)::Union{Missing, Float64}
    n = Float64(num)
    d = Float64(den)
    if !isfinite(n) || !isfinite(d) || d <= 0.0
        return missing
    end
    return n / d
end

@inline function _safe_ratio(num::Missing, den)::Union{Missing, Float64}
    return missing
end

@inline function _safe_ratio(num, den::Missing)::Union{Missing, Float64}
    return missing
end

@inline function _safe_ratio(num, den)::Union{Missing, Float64}
    if !(num isa Real) || !(den isa Real)
        return missing
    end
    return _safe_ratio(num::Real, den::Real)
end

@inline function _key_token(value)::String
    if value isa Missing || value === nothing
        return "_"
    elseif value isa AbstractFloat
        v = Float64(value)
        return isfinite(v) ? string(round(v; digits=6)) : "_"
    end
    return string(value)
end

@inline function _mission_family(row)::String
    category = hasproperty(row, :category) ? lowercase(string(getproperty(row, :category))) : ""
    scenario = hasproperty(row, :scenario) ? lowercase(string(getproperty(row, :scenario))) : ""
    mission_time_s = hasproperty(row, :mission_time_s) && getproperty(row, :mission_time_s) isa Real ? Float64(getproperty(row, :mission_time_s)) : 0.0
    satellites = hasproperty(row, :satellites) && getproperty(row, :satellites) isa Integer ? Int(getproperty(row, :satellites)) : 1
    dynamic_effectors = hasproperty(row, :dynamic_effectors) ? lowercase(string(getproperty(row, :dynamic_effectors))) : ""

    if category == "montecarlo" || occursin("montecarlo", scenario)
        return mission_time_s <= 7200.0 ? "mc_short" : "mc_long"
    end
    if satellites >= 2
        return "multi_sat"
    end
    if occursin("nbody", dynamic_effectors) || occursin("harmonic", dynamic_effectors)
        return "high_fidelity_nbody"
    end
    if mission_time_s > 7200.0 || category == "mission_length"
        return "long"
    end
    return "short_light"
end

@inline function _match_key(row)::String
    parts = (
        _key_token(hasproperty(row, :pass) ? getproperty(row, :pass) : missing),
        _key_token(hasproperty(row, :scenario) ? getproperty(row, :scenario) : missing),
        _key_token(hasproperty(row, :category) ? getproperty(row, :category) : missing),
        _key_token(hasproperty(row, :repeat) ? getproperty(row, :repeat) : missing),
        _key_token(hasproperty(row, :seed) ? getproperty(row, :seed) : missing),
        _key_token(hasproperty(row, :mission_time_s) ? getproperty(row, :mission_time_s) : missing),
        _key_token(hasproperty(row, :satellites) ? getproperty(row, :satellites) : missing)
    )
    return join(parts, "|")
end

function _baseline_artifact(artifacts::Vector{ModeRunArtifacts})::ModeRunArtifacts
    idx = findfirst(a -> a.mode == :serial, artifacts)
    idx === nothing && error("Smart ladder requires a serial baseline rung with mode=:serial.")
    return artifacts[idx]
end

function _prepare_speed_sample_table(raw_df::DataFrame)::DataFrame
    rows = NamedTuple[]
    for row in eachrow(raw_df)
        if hasproperty(row, :solve_success) && row.solve_success !== true
            continue
        end
        total_time = hasproperty(row, :total_time_s) ? getproperty(row, :total_time_s) : missing
        if !(total_time isa Real) || !isfinite(Float64(total_time)) || Float64(total_time) <= 0.0
            continue
        end
        push!(rows, (
            match_key=_match_key(row),
            total_time_s=Float64(total_time),
            mission_family=_mission_family(row)
        ))
    end
    sample_df = DataFrame(rows)
    nrow(sample_df) == 0 && return sample_df
    return combine(
        groupby(sample_df, :match_key),
        :total_time_s => mean => :total_time_s,
        :mission_family => (v -> first(v)) => :mission_family
    )
end

function _build_vs_r0_speedup_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact(artifacts)
    rows = NamedTuple[]
    for artifact in artifacts
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            backend=artifact.backend,
            total_elapsed_s=artifact.elapsed_s,
            run_benchmarks_elapsed_s=artifact.bench_elapsed_s,
            run_per_orbit_elapsed_s=artifact.orbit_elapsed_s,
            total_speedup_vs_r0=_safe_ratio(baseline.elapsed_s, artifact.elapsed_s),
            run_benchmarks_speedup_vs_r0=_safe_ratio(baseline.bench_elapsed_s, artifact.bench_elapsed_s),
            run_per_orbit_speedup_vs_r0=_safe_ratio(baseline.orbit_elapsed_s, artifact.orbit_elapsed_s)
        ))
    end
    df = DataFrame(rows)
    sort!(df, :mode)
    return df
end

function _build_mission_family_speedup_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact(artifacts)
    baseline_samples = _prepare_speed_sample_table(baseline.raw_df)
    if nrow(baseline_samples) == 0
        return DataFrame()
    end
    baseline_samples = select(baseline_samples, :match_key, :total_time_s => :r0_total_time_s)

    rows = NamedTuple[]
    for artifact in artifacts
        artifact.mode == :serial && continue
        rung_samples = _prepare_speed_sample_table(artifact.raw_df)
        nrow(rung_samples) == 0 && continue
        joined = innerjoin(
            baseline_samples,
            select(rung_samples, :match_key, :mission_family, :total_time_s => :rung_total_time_s),
            on=:match_key
        )
        nrow(joined) == 0 && continue
        joined[!, :speedup_vs_r0] = [Float64(r0) / Float64(rt) for (r0, rt) in zip(joined.r0_total_time_s, joined.rung_total_time_s)]
        joined = joined[isfinite.(joined.speedup_vs_r0) .& (joined.speedup_vs_r0 .> 0.0), :]
        nrow(joined) == 0 && continue

        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        for grp in groupby(joined, :mission_family)
            vals = Float64.(grp.speedup_vs_r0)
            sample_count = length(vals)
            slower_count = count(v -> v < 1.0, vals)
            push!(rows, (
                rung=rung_label,
                mode=String(artifact.mode),
                mission_family=String(first(grp.mission_family)),
                samples=sample_count,
                median_speedup_vs_r0=median(vals),
                p90_speedup_vs_r0=quantile(vals, 0.9),
                worst_slowdown_x=maximum(1.0 ./ vals),
                slower_share_pct=100.0 * slower_count / sample_count
            ))
        end
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, [:rung, :mission_family])
    end
    return df
end

@inline function _mean_skipmissing(values)
    vec = collect(skipmissing(values))
    isempty(vec) && return missing
    return mean(Float64.(vec))
end

@inline function _safe_rel_error(r0, value)
    if !(r0 isa Real) || !(value isa Real)
        return missing
    end
    r0_f = Float64(r0)
    value_f = Float64(value)
    if !isfinite(r0_f) || !isfinite(value_f)
        return missing
    end
    denom = max(abs(r0_f), eps(Float64))
    return abs(value_f - r0_f) / denom
end

@inline function _error_stats(values)::NamedTuple
    vec = collect(skipmissing(values))
    isempty(vec) && return (median=missing, p90=missing, max=missing)
    return (
        median=median(vec),
        p90=quantile(vec, 0.9),
        max=maximum(vec)
    )
end

function _prepare_accuracy_sample_table(raw_df::DataFrame)::DataFrame
    required = (
        :final_primary_pos_norm_m,
        :final_primary_vel_norm_mps,
        :final_primary_mass_kg
    )
    for column in required
        _has_column(raw_df, column) || return DataFrame()
    end

    rows = NamedTuple[]
    for row in eachrow(raw_df)
        if hasproperty(row, :solve_success) && row.solve_success !== true
            continue
        end
        push!(rows, (
            match_key=_match_key(row),
            final_primary_pos_norm_m=getproperty(row, :final_primary_pos_norm_m),
            final_primary_vel_norm_mps=getproperty(row, :final_primary_vel_norm_mps),
            final_primary_mass_kg=getproperty(row, :final_primary_mass_kg)
        ))
    end

    sample_df = DataFrame(rows)
    nrow(sample_df) == 0 && return sample_df
    return combine(
        groupby(sample_df, :match_key),
        :final_primary_pos_norm_m => _mean_skipmissing => :final_primary_pos_norm_m,
        :final_primary_vel_norm_mps => _mean_skipmissing => :final_primary_vel_norm_mps,
        :final_primary_mass_kg => _mean_skipmissing => :final_primary_mass_kg
    )
end

function _build_accuracy_parity_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    baseline = _baseline_artifact(artifacts)
    baseline_samples = _prepare_accuracy_sample_table(baseline.raw_df)
    if nrow(baseline_samples) == 0
        return DataFrame()
    end
    baseline_samples = select(
        baseline_samples,
        :match_key,
        :final_primary_pos_norm_m => :r0_final_primary_pos_norm_m,
        :final_primary_vel_norm_mps => :r0_final_primary_vel_norm_mps,
        :final_primary_mass_kg => :r0_final_primary_mass_kg
    )

    rows = NamedTuple[]
    for artifact in artifacts
        rung_samples = _prepare_accuracy_sample_table(artifact.raw_df)
        nrow(rung_samples) == 0 && continue
        joined = innerjoin(
            baseline_samples,
            select(
                rung_samples,
                :match_key,
                :final_primary_pos_norm_m => :rung_final_primary_pos_norm_m,
                :final_primary_vel_norm_mps => :rung_final_primary_vel_norm_mps,
                :final_primary_mass_kg => :rung_final_primary_mass_kg
            ),
            on=:match_key
        )
        nrow(joined) == 0 && continue

        pos_err = [_safe_rel_error(r0, rt) for (r0, rt) in zip(joined.r0_final_primary_pos_norm_m, joined.rung_final_primary_pos_norm_m)]
        vel_err = [_safe_rel_error(r0, rt) for (r0, rt) in zip(joined.r0_final_primary_vel_norm_mps, joined.rung_final_primary_vel_norm_mps)]
        mass_err = [_safe_rel_error(r0, rt) for (r0, rt) in zip(joined.r0_final_primary_mass_kg, joined.rung_final_primary_mass_kg)]

        pos_stats = _error_stats(pos_err)
        vel_stats = _error_stats(vel_err)
        mass_stats = _error_stats(mass_err)
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            samples=nrow(joined),
            pos_rel_err_median_pct=ismissing(pos_stats.median) ? missing : 100.0 * pos_stats.median,
            pos_rel_err_p90_pct=ismissing(pos_stats.p90) ? missing : 100.0 * pos_stats.p90,
            pos_rel_err_max_pct=ismissing(pos_stats.max) ? missing : 100.0 * pos_stats.max,
            vel_rel_err_median_pct=ismissing(vel_stats.median) ? missing : 100.0 * vel_stats.median,
            vel_rel_err_p90_pct=ismissing(vel_stats.p90) ? missing : 100.0 * vel_stats.p90,
            vel_rel_err_max_pct=ismissing(vel_stats.max) ? missing : 100.0 * vel_stats.max,
            mass_rel_err_median_pct=ismissing(mass_stats.median) ? missing : 100.0 * mass_stats.median,
            mass_rel_err_p90_pct=ismissing(mass_stats.p90) ? missing : 100.0 * mass_stats.p90,
            mass_rel_err_max_pct=ismissing(mass_stats.max) ? missing : 100.0 * mass_stats.max
        ))
    end

    df = DataFrame(rows)
    if nrow(df) > 0
        sort!(df, :mode)
    end
    return df
end

function _build_route_mix_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String}
)::DataFrame
    rows = NamedTuple[]
    for artifact in artifacts
        raw_df = artifact.raw_df
        total_rows = nrow(raw_df)
        successful_rows = _has_column(raw_df, :solve_success) ? count(==(true), raw_df.solve_success) : total_rows
        route_vals = _has_column(raw_df, :outer_route) ? lowercase.(string.(raw_df.outer_route)) : fill("unknown", total_rows)
        none_count = count(==("none"), route_vals)
        threads_count = count(==("threads"), route_vals)
        process_count = count(==("process"), route_vals)
        other_count = max(0, total_rows - (none_count + threads_count + process_count))
        denom = total_rows > 0 ? total_rows : 1
        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            total_rows=total_rows,
            successful_rows=successful_rows,
            none_count=none_count,
            threads_count=threads_count,
            process_count=process_count,
            other_count=other_count,
            none_pct=100.0 * none_count / denom,
            threads_pct=100.0 * threads_count / denom,
            process_pct=100.0 * process_count / denom,
            other_pct=100.0 * other_count / denom
        ))
    end
    df = DataFrame(rows)
    sort!(df, :mode)
    return df
end

function _gate_counts(df::Union{Nothing, DataFrame})::NamedTuple
    if df === nothing
        return (rows=0, pass_rows=0, pass_rate_pct=missing)
    end
    rows = nrow(df)
    if rows == 0 || !_has_column(df, :pass_all)
        return (rows=rows, pass_rows=0, pass_rate_pct=missing)
    end
    pass_rows = count(Bool.(df.pass_all))
    return (rows=rows, pass_rows=pass_rows, pass_rate_pct=100.0 * pass_rows / rows)
end

function _build_fidelity_parity_table(
    artifacts::Vector{ModeRunArtifacts},
    rung_label_by_mode::Dict{Symbol, String},
    config::SmartLadderConfig
)::DataFrame
    rows = NamedTuple[]
    for artifact in artifacts
        raw_df = artifact.raw_df
        total_rows = nrow(raw_df)
        success_rows = _has_column(raw_df, :solve_success) ? count(==(true), raw_df.solve_success) : total_rows
        success_rate = total_rows > 0 ? 100.0 * success_rows / total_rows : missing
        split_counts = _gate_counts(artifact.split_gate_df)

        agg_dir = dirname(artifact.raw_path)
        multirate_csv_path = _latest_artifact_path_optional(
            agg_dir,
            "multirate_rollout_gate_agg",
            config.profile.name,
            ".csv"
        )
        multirate_df = isnothing(multirate_csv_path) ? nothing : CSV.read(multirate_csv_path, DataFrame)
        multirate_counts = _gate_counts(multirate_df)

        rung_label = get(rung_label_by_mode, artifact.mode, string(artifact.mode))
        push!(rows, (
            rung=rung_label,
            mode=String(artifact.mode),
            sample_rows=total_rows,
            success_rows=success_rows,
            success_rate_pct=success_rate,
            split_gate_rows=split_counts.rows,
            split_gate_pass_rows=split_counts.pass_rows,
            split_gate_pass_rate_pct=split_counts.pass_rate_pct,
            multirate_gate_rows=multirate_counts.rows,
            multirate_gate_pass_rows=multirate_counts.pass_rows,
            multirate_gate_pass_rate_pct=multirate_counts.pass_rate_pct
        ))
    end
    df = DataFrame(rows)
    sort!(df, :mode)
    return df
end

function _write_aggregate_rung_report(
    path::String,
    config::SmartLadderConfig,
    rung::LadderRungSpec,
    runs::Vector{LadderPassResult},
    bench_elapsed_s::Float64,
    orbit_elapsed_s::Float64,
    total_elapsed_s::Float64;
    split_gate_csv_path::Union{Nothing, String}=nothing,
    split_gate_pass_rows::Int=0,
    split_gate_total_rows::Int=0,
    multirate_gate_csv_path::Union{Nothing, String}=nothing,
    multirate_gate_pass_rows::Int=0,
    multirate_gate_total_rows::Int=0
)
    open(path, "w") do io
        println(io, "# Smart Parallel Ladder Aggregated Rung Report")
        println(io)
        println(io, "- Generated (UTC): $(now(UTC))")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Rung: `$(rung.label)` (`$(rung.mode)`) ")
        println(io, "- Description: $(rung.description)")
        println(io, "- Matrix: `$(rung.matrix)`")
        println(io, "- Backend: `$(rung.backend)`")
        println(io, "- Inner adaptive: `$(rung.inner_adaptive)`")
        println(io, "- Outer-route adaptive: `$(rung.outer_route_adaptive)`")
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
        if multirate_gate_total_rows > 0
            println(io, "- Multirate rollout gate pass rows: `$(multirate_gate_pass_rows)/$(multirate_gate_total_rows)`")
        else
            println(io, "- Multirate rollout gate rows: `0`")
        end
        if !(multirate_gate_csv_path === nothing)
            println(io, "- Multirate rollout gate aggregated CSV: `$(multirate_gate_csv_path)`")
        end
        println(io)
        println(io, "## Source Runs")
        println(io)
        for run in runs
            println(io, "- pass=$(run.pass): raw=`$(run.artifact.raw_path)`, per-orbit raw=`$(run.artifact.orbit_raw_path)`, report=`$(run.artifact.report_path)`")
        end
    end
end

function _aggregate_rung_artifacts(
    config::SmartLadderConfig,
    rung::LadderRungSpec,
    runs::Vector{LadderPassResult}
)::ModeRunArtifacts
    isempty(runs) && error("No pass data collected for rung '$(rung.label)'.")

    raw_df = DataFrame()
    orbit_raw_df = DataFrame()
    split_gate_df = DataFrame()
    multirate_gate_df = DataFrame()
    for run in runs
        raw_df = vcat(raw_df, _tag_rung_column(run.artifact.raw_df, rung; pass_idx=run.pass); cols=:union)
        local_orbit_raw = CSV.read(run.artifact.orbit_raw_path, DataFrame)
        orbit_raw_df = vcat(orbit_raw_df, _tag_rung_column(local_orbit_raw, rung; pass_idx=run.pass); cols=:union)
        if !(run.artifact.split_gate_df === nothing)
            split_gate_df = vcat(split_gate_df, _tag_rung_column(run.artifact.split_gate_df, rung; pass_idx=run.pass); cols=:union)
        end
        run_outdir = dirname(run.artifact.raw_path)
        run_multirate_csv_path = _latest_artifact_path_optional(run_outdir, "multirate_rollout_gate", config.profile.name, ".csv")
        if !(run_multirate_csv_path === nothing)
            run_multirate_df = CSV.read(run_multirate_csv_path, DataFrame)
            multirate_gate_df = vcat(multirate_gate_df, _tag_rung_column(run_multirate_df, rung; pass_idx=run.pass); cols=:union)
        end
    end

    summary_df = summarize_results(raw_df)
    orbit_summary_df = summarize_per_orbit_results(orbit_raw_df)

    bench_elapsed_s = mean([run.artifact.bench_elapsed_s for run in runs])
    split_gate_elapsed_s = mean([run.artifact.split_gate_elapsed_s for run in runs])
    orbit_elapsed_s = mean([run.artifact.orbit_elapsed_s for run in runs])
    total_elapsed_s = bench_elapsed_s + split_gate_elapsed_s + orbit_elapsed_s

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    agg_outdir = joinpath(config.outdir, "aggregate", rung.label)
    mkpath(agg_outdir)

    raw_path = joinpath(agg_outdir, "runtime_raw_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    summary_path = joinpath(agg_outdir, "runtime_summary_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    orbit_raw_path = joinpath(agg_outdir, "runtime_per_orbit_raw_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    orbit_summary_path = joinpath(agg_outdir, "runtime_per_orbit_summary_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    stage_timing_path = joinpath(agg_outdir, "runtime_stage_timing_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    hardware_info_path = joinpath(agg_outdir, "runtime_hardware_info_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv")
    split_gate_csv_path = nrow(split_gate_df) > 0 ? joinpath(agg_outdir, "split_rollout_gate_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv") : nothing
    multirate_gate_csv_path = nrow(multirate_gate_df) > 0 ? joinpath(agg_outdir, "multirate_rollout_gate_agg_$(config.profile.name)_$(rung.label)_$(stamp).csv") : nothing
    report_path = joinpath(agg_outdir, "runtime_report_agg_$(config.profile.name)_$(rung.label)_$(stamp).md")

    split_gate_pass_rows = (nrow(split_gate_df) > 0 && _has_column(split_gate_df, :pass_all)) ? count(Bool.(split_gate_df.pass_all)) : 0
    split_gate_total_rows = nrow(split_gate_df)
    multirate_gate_pass_rows = (nrow(multirate_gate_df) > 0 && _has_column(multirate_gate_df, :pass_all)) ? count(Bool.(multirate_gate_df.pass_all)) : 0
    multirate_gate_total_rows = nrow(multirate_gate_df)

    hw = _runtime_hardware_snapshot()
    hardware_info_df = DataFrame([
        (
            profile=config.profile.name,
            rung=rung.label,
            rung_mode=String(rung.mode),
            matrix=String(rung.matrix),
            backend=rung.backend,
            inner_adaptive=rung.inner_adaptive,
            outer_route_adaptive=rung.outer_route_adaptive,
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
    if !(multirate_gate_csv_path === nothing)
        CSV.write(multirate_gate_csv_path, multirate_gate_df)
    end

    _write_aggregate_rung_report(
        report_path,
        config,
        rung,
        runs,
        bench_elapsed_s,
        orbit_elapsed_s,
        total_elapsed_s;
        split_gate_csv_path=split_gate_csv_path,
        split_gate_pass_rows=split_gate_pass_rows,
        split_gate_total_rows=split_gate_total_rows,
        multirate_gate_csv_path=multirate_gate_csv_path,
        multirate_gate_pass_rows=multirate_gate_pass_rows,
        multirate_gate_total_rows=multirate_gate_total_rows
    )

    backend_label =
        "$(rung.backend)|matrix=$(rung.matrix)|inner=$(rung.inner_adaptive ? "adaptive" : "static")|outer=$(rung.outer_route_adaptive ? "adaptive" : "static")"
    return ModeRunArtifacts(
        mode=rung.mode,
        backend=backend_label,
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

function _write_smart_ladder_report(
    path::String,
    config::SmartLadderConfig,
    rungs::Vector{LadderRungSpec},
    overview_df::DataFrame,
    comparison_df::DataFrame,
    speedup_vs_r0_df::DataFrame,
    mission_family_df::DataFrame,
    fidelity_df::DataFrame,
    accuracy_df::DataFrame,
    route_mix_df::DataFrame,
    artifacts::Vector{ModeRunArtifacts},
    run_order_df::DataFrame;
    comparison_plot_paths::Vector{String}=String[],
    comparison_csv_path::String="",
    overview_csv_path::String="",
    speedup_vs_r0_csv_path::String="",
    mission_family_csv_path::String="",
    fidelity_csv_path::String="",
    accuracy_csv_path::String="",
    route_mix_csv_path::String="",
    merged_raw_csv_path::String="",
    merged_summary_csv_path::String="",
    run_order_csv_path::String=""
)
    generated = string(now(UTC))
    nthreads = Threads.nthreads()
    cpu_threads = Sys.CPU_THREADS
    rung_labels = join([r.label for r in rungs], ", ")

    open(path, "w") do io
        println(io, "# SpaceAGORA Smart Parallel Ladder")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- Threads in wrapper process: `$nthreads`")
        println(io, "- Detected CPU threads: `$cpu_threads`")
        println(io, "- Profile: `$(config.profile.name)`")
        println(io, "- Rungs executed: `$rung_labels`")
        println(io, "- Passes: `$(config.passes)`")
        println(io, "- Randomized rung order: `$(config.randomize_rung_order)` (seed=`$(config.random_seed)`)")
        println(io, "- Outer-only backend: `$(config.outer_only_backend)`")
        println(io, "- Include control_stress in per-orbit stage: `$(config.include_control_stress_per_orbit)`")
        println(io, "- Control_stress full schedule override: warmup=`$(config.control_stress_warmup_full)`, repeats=`$(config.control_stress_repeats_full)`")
        println(io)

        println(io, "## Ladder Definition")
        println(io)
        definition_df = DataFrame(
            rung=[r.label for r in rungs],
            mode=[String(r.mode) for r in rungs],
            matrix=[String(r.matrix) for r in rungs],
            backend=[r.backend for r in rungs],
            inner_adaptive=[r.inner_adaptive for r in rungs],
            outer_route_adaptive=[r.outer_route_adaptive for r in rungs],
            description=[r.description for r in rungs]
        )
        _write_markdown_table(io, definition_df)
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
        println(io, "_Note: `serial_total_time_mean_s` corresponds to rung `r0_true_serial`._")
        println(io)
        _write_markdown_table(io, comparison_df)
        println(io)

        println(io, "## Speedup vs R0 (Wall Time and Stages)")
        println(io)
        _write_markdown_table(io, speedup_vs_r0_df)
        println(io)

        println(io, "## Mission-Family Speedup Distribution vs R0")
        println(io)
        println(io, "_Family rules: `mc_short`/`mc_long` by mission time <=/> 7200 s; `multi_sat` for satellites >= 2; `high_fidelity_nbody` for N-body/harmonics; `long` for mission_time > 7200 s; otherwise `short_light`._")
        println(io)
        if nrow(mission_family_df) > 0
            _write_markdown_table(io, mission_family_df)
        else
            println(io, "- No mission-family speedup samples available.")
        end
        println(io)

        println(io, "## Fidelity Parity")
        println(io)
        _write_markdown_table(io, fidelity_df)
        println(io)

        println(io, "## Accuracy Parity vs R0")
        println(io)
        println(io, "_Metrics are relative errors (%) on terminal primary-spacecraft norms (position, velocity, mass) against `r0_true_serial`._")
        println(io)
        if nrow(accuracy_df) > 0
            _write_markdown_table(io, accuracy_df)
        else
            println(io, "- No accuracy parity samples available.")
        end
        println(io)

        println(io, "## Route Mix")
        println(io)
        _write_markdown_table(io, route_mix_df)
        println(io)

        println(io, "## Plots")
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
        println(io, "- Speedup vs R0 CSV: `$(speedup_vs_r0_csv_path)`")
        println(io, "- Mission-family speedup CSV: `$(mission_family_csv_path)`")
        println(io, "- Fidelity parity CSV: `$(fidelity_csv_path)`")
        println(io, "- Accuracy parity CSV: `$(accuracy_csv_path)`")
        println(io, "- Route mix CSV: `$(route_mix_csv_path)`")
        println(io, "- Merged raw CSV (all passes): `$(merged_raw_csv_path)`")
        println(io, "- Merged summary CSV (all passes): `$(merged_summary_csv_path)`")
        println(io, "- Run order CSV: `$(run_order_csv_path)`")
        for artifact in artifacts
            println(io, "- Rung `$(artifact.mode)` aggregated raw: `$(artifact.raw_path)`")
            println(io, "- Rung `$(artifact.mode)` aggregated summary: `$(artifact.summary_path)`")
            println(io, "- Rung `$(artifact.mode)` aggregated per-orbit summary: `$(artifact.orbit_summary_path)`")
            if !isempty(artifact.stage_timing_path)
                println(io, "- Rung `$(artifact.mode)` aggregated stage timing: `$(artifact.stage_timing_path)`")
            end
            if !isempty(artifact.hardware_info_path)
                println(io, "- Rung `$(artifact.mode)` aggregated hardware info: `$(artifact.hardware_info_path)`")
            end
            if !(artifact.split_gate_csv_path === nothing)
                println(io, "- Rung `$(artifact.mode)` aggregated split gate CSV: `$(artifact.split_gate_csv_path)`")
            end
            println(io, "- Rung `$(artifact.mode)` aggregated report: `$(artifact.report_path)`")
        end
        println(io)

        println(io, "## Reproducibility")
        println(io, "```bash")
        process_suffix = isnothing(config.process_workers) ? "" : " --process-workers=$(config.process_workers)"
        println(
            io,
            "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=.AGORA test/performance_smart_parallel_ladder.jl " *
            "--profile=$(config.profile.name) --outdir=$(config.outdir) --clean=1 --passes=$(config.passes) " *
            "--randomize-rung-order=$(config.randomize_rung_order ? 1 : 0) --seed=$(config.random_seed) " *
            "--outer-only-backend=$(config.outer_only_backend)$(process_suffix) " *
            "--include-control-stress-per-orbit=$(config.include_control_stress_per_orbit ? 1 : 0) " *
            "--control-stress-repeats-full=$(config.control_stress_repeats_full) " *
            "--control-stress-warmup-full=$(config.control_stress_warmup_full)"
        )
        println(io, "```")
    end
end

function main_smart_parallel_ladder()
    config = parse_smart_ladder_cli()
    if config.clean
        rm(config.outdir; recursive=true, force=true)
    end
    mkpath(config.outdir)

    rungs = _ladder_rungs(config)
    rng = MersenneTwister(config.random_seed)

    println("Smart parallel ladder profile: $(config.profile.name)")
    println("Outdir: $(config.outdir)")
    println("Clean outdir: $(config.clean)")
    println("Rungs: $(join([r.label for r in rungs], ", "))")
    println("Passes: $(config.passes)")
    println("Randomized rung order: $(config.randomize_rung_order) (seed=$(config.random_seed))")
    println("Outer-only backend: $(config.outer_only_backend)")
    println("Include control_stress in per-orbit stage: $(config.include_control_stress_per_orbit)")
    println("Control_stress full schedule override: warmup=$(config.control_stress_warmup_full), repeats=$(config.control_stress_repeats_full)")
    println("Wrapper Threads.nthreads()=$(Threads.nthreads()), Sys.CPU_THREADS=$(Sys.CPU_THREADS)")

    pass_results = LadderPassResult[]
    order_rows = NamedTuple[]
    for pass_idx in 1:config.passes
        pass_outdir = joinpath(config.outdir, "pass_$(lpad(pass_idx, 2, '0'))")
        mkpath(pass_outdir)
        ordered = _ordered_rungs(rungs, rng, config.randomize_rung_order)
        order_line = join([r.label for r in ordered], " -> ")
        println("[smart-ladder] pass=$(pass_idx) rung-order=$(order_line)")
        for (order_idx, rung) in enumerate(ordered)
            push!(order_rows, (
                pass=pass_idx,
                order_index=order_idx,
                rung=rung.label,
                mode=String(rung.mode),
                matrix=String(rung.matrix),
                backend=rung.backend,
                inner_adaptive=rung.inner_adaptive,
                outer_route_adaptive=rung.outer_route_adaptive
            ))
            artifact = run_rung(rung, config, pass_idx, pass_outdir)
            push!(pass_results, LadderPassResult(pass=pass_idx, rung=rung, artifact=artifact))
        end
    end

    aggregated_artifacts = ModeRunArtifacts[]
    for rung in rungs
        runs = [result for result in pass_results if result.rung.label == rung.label]
        push!(aggregated_artifacts, _aggregate_rung_artifacts(config, rung, runs))
    end

    comparison_df = build_comparison_table(aggregated_artifacts)
    overview_df = build_mode_overview(aggregated_artifacts)
    rung_map = Dict(r.mode => r for r in rungs)
    rung_label_by_mode = Dict(r.mode => r.label for r in rungs)
    overview_df.mode = [
        haskey(rung_map, Symbol(mode_name)) ? rung_map[Symbol(mode_name)].label : String(mode_name)
        for mode_name in String.(overview_df.mode)
    ]
    speedup_vs_r0_df = _build_vs_r0_speedup_table(aggregated_artifacts, rung_label_by_mode)
    mission_family_df = _build_mission_family_speedup_table(aggregated_artifacts, rung_label_by_mode)
    fidelity_df = _build_fidelity_parity_table(aggregated_artifacts, rung_label_by_mode, config)
    accuracy_df = _build_accuracy_parity_table(aggregated_artifacts, rung_label_by_mode)
    route_mix_df = _build_route_mix_table(aggregated_artifacts, rung_label_by_mode)

    merged_raw_df = DataFrame()
    merged_summary_df = DataFrame()
    for result in pass_results
        merged_raw_df = vcat(
            merged_raw_df,
            _tag_rung_column(result.artifact.raw_df, result.rung; pass_idx=result.pass);
            cols=:union
        )
        merged_summary_df = vcat(
            merged_summary_df,
            _tag_rung_column(result.artifact.summary_df, result.rung; pass_idx=result.pass);
            cols=:union
        )
    end

    run_order_df = DataFrame(order_rows)

    stamp = Dates.format(now(UTC), dateformat"yyyymmdd_HHMMSS")
    comparison_path = joinpath(config.outdir, "smart_parallel_ladder_compare_summary_$(config.profile.name)_$(stamp).csv")
    overview_path = joinpath(config.outdir, "smart_parallel_ladder_mode_overview_$(config.profile.name)_$(stamp).csv")
    speedup_vs_r0_path = joinpath(config.outdir, "smart_parallel_ladder_speedup_vs_r0_$(config.profile.name)_$(stamp).csv")
    mission_family_path = joinpath(config.outdir, "smart_parallel_ladder_mission_family_speedup_$(config.profile.name)_$(stamp).csv")
    fidelity_path = joinpath(config.outdir, "smart_parallel_ladder_fidelity_parity_$(config.profile.name)_$(stamp).csv")
    accuracy_path = joinpath(config.outdir, "smart_parallel_ladder_accuracy_parity_$(config.profile.name)_$(stamp).csv")
    route_mix_path = joinpath(config.outdir, "smart_parallel_ladder_route_mix_$(config.profile.name)_$(stamp).csv")
    merged_raw_path = joinpath(config.outdir, "smart_parallel_ladder_raw_merged_$(config.profile.name)_$(stamp).csv")
    merged_summary_path = joinpath(config.outdir, "smart_parallel_ladder_summary_merged_$(config.profile.name)_$(stamp).csv")
    run_order_path = joinpath(config.outdir, "smart_parallel_ladder_run_order_$(config.profile.name)_$(stamp).csv")
    report_path = joinpath(config.outdir, "smart_parallel_ladder_report_$(config.profile.name)_$(stamp).md")
    comparison_plot_paths = generate_pipeline_comparison_plots(
        config.outdir,
        config.profile,
        "smart_parallel_ladder_$(stamp)",
        overview_df,
        comparison_df
    )

    CSV.write(comparison_path, comparison_df)
    CSV.write(overview_path, overview_df)
    CSV.write(speedup_vs_r0_path, speedup_vs_r0_df)
    CSV.write(mission_family_path, mission_family_df)
    CSV.write(fidelity_path, fidelity_df)
    CSV.write(accuracy_path, accuracy_df)
    CSV.write(route_mix_path, route_mix_df)
    CSV.write(merged_raw_path, merged_raw_df)
    CSV.write(merged_summary_path, merged_summary_df)
    CSV.write(run_order_path, run_order_df)
    _write_smart_ladder_report(
        report_path,
        config,
        rungs,
        overview_df,
        comparison_df,
        speedup_vs_r0_df,
        mission_family_df,
        fidelity_df,
        accuracy_df,
        route_mix_df,
        aggregated_artifacts,
        run_order_df;
        comparison_plot_paths=comparison_plot_paths,
        comparison_csv_path=comparison_path,
        overview_csv_path=overview_path,
        speedup_vs_r0_csv_path=speedup_vs_r0_path,
        mission_family_csv_path=mission_family_path,
        fidelity_csv_path=fidelity_path,
        accuracy_csv_path=accuracy_path,
        route_mix_csv_path=route_mix_path,
        merged_raw_csv_path=merged_raw_path,
        merged_summary_csv_path=merged_summary_path,
        run_order_csv_path=run_order_path
    )

    println()
    println("Smart parallel ladder complete.")
    println("comparison summary: $(comparison_path)")
    println("mode overview: $(overview_path)")
    println("speedup vs r0: $(speedup_vs_r0_path)")
    println("mission-family speedup: $(mission_family_path)")
    println("fidelity parity: $(fidelity_path)")
    println("accuracy parity: $(accuracy_path)")
    println("route mix: $(route_mix_path)")
    println("merged raw: $(merged_raw_path)")
    println("merged summary: $(merged_summary_path)")
    println("run order: $(run_order_path)")
    println("comparison plots: $(length(comparison_plot_paths))")
    println("report: $(report_path)")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_smart_parallel_ladder()
end
