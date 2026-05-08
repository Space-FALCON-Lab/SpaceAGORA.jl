const _SMART_LADDER_STUDIES_DIR = dirname(@__DIR__)
const _SMART_LADDER_BENCHMARKS_DIR = dirname(_SMART_LADDER_STUDIES_DIR)
include(joinpath(_SMART_LADDER_BENCHMARKS_DIR, "scripts", "performance_paper_pipeline.jl"))

using Random
using SHA

const SMART_LADDER_DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "performance", "smart_parallel_ladder")
const SMART_LADDER_RUNTIME_SCRIPT = joinpath(_SMART_LADDER_STUDIES_DIR, "performance_runtime_analysis.jl")

@inline function _smart_ladder_smoke_mode()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_SMART_LADDER_SMOKE", "0")))
    return raw in ("1", "true", "yes", "on")
end

@inline function _is_julia_project_dir(path::String)::Bool
    return isdir(path) && isfile(joinpath(path, "Project.toml"))
end

function _resolve_smart_ladder_project()::String
    candidates = Pair{String, String}[]
    env_project = strip(get(ENV, "SPACEAGORA_SMART_LADDER_PROJECT", ""))
    if !isempty(env_project)
        push!(candidates, "SPACEAGORA_SMART_LADDER_PROJECT" => env_project)
    end
    push!(candidates, "repo_dot_AGORA" => joinpath(REPO_ROOT, ".AGORA"))
    push!(candidates, "repo_root" => REPO_ROOT)

    checked = String[]
    for (label, raw_path) in candidates
        path = abspath(raw_path)
        push!(checked, "$(label)=$(path)")
        if _is_julia_project_dir(path)
            return path
        end
    end

    error(
        "Unable to resolve smart ladder project path. Checked: $(join(checked, ", ")). " *
        "Set SPACEAGORA_SMART_LADDER_PROJECT to a valid Julia project directory (must contain Project.toml)."
    )
end

const SMART_LADDER_PROJECT = _resolve_smart_ladder_project()

Base.@kwdef struct SmartLadderConfig
    profile::ProfileSpec
    outdir::String
    clean::Bool
    passes::Int
    randomize_rung_order::Bool
    random_seed::Int
    outer_only_backend::String
    process_workers::Union{Nothing, Int}
    include_layer_attribution::Bool
    rung_filter::Vector{String}
    include_control_stress_per_orbit::Bool
    control_stress_repeats_full::Int
    control_stress_warmup_full::Int
    solver_axis::Symbol
    solver_mode::String
    solver_factor_modes::Vector{String}
    trajectory_output::Bool
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
    solver_label::String
    solver_mode::Union{Nothing, String}
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

@inline function _parse_solver_axis(raw::AbstractString)::Symbol
    token = lowercase(strip(String(raw)))
    token in ("inherit", "frozen", "factorial") || throw(
        ArgumentError("solver-axis must be one of: inherit, frozen, factorial (got '$raw').")
    )
    return Symbol(token)
end

@inline function _parse_solver_factor_modes(raw::AbstractString)::Vector{String}
    token = strip(String(raw))
    isempty(token) && return String[]
    values = String[]
    seen = Set{String}()
    for part in split(token, ",")
        mode = lowercase(strip(part))
        isempty(mode) && continue
        if !(mode in seen)
            push!(values, mode)
            push!(seen, mode)
        end
    end
    return values
end

function _parse_rung_filter(raw::AbstractString)::Vector{String}
    token = strip(String(raw))
    isempty(token) && return String[]
    values = String[]
    seen = Set{String}()
    for part in split(token, ",")
        value = strip(part)
        isempty(value) && continue
        if !(value in seen)
            push!(values, value)
            push!(seen, value)
        end
    end
    return values
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
    include_layer_attribution = _parse_bool_token(get(ENV, "SPACEAGORA_SMART_LADDER_LAYER_ATTRIBUTION", "1"))
    rung_filter = _parse_rung_filter(get(ENV, "SPACEAGORA_SMART_LADDER_RUNGS", ""))
    include_control_stress_per_orbit = _parse_bool_token(get(ENV, "SPACEAGORA_PERF_INCLUDE_CONTROL_STRESS_PER_ORBIT", "1"))
    control_stress_repeats_full = _parse_positive_int(
        get(ENV, "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL", "3"),
        "SPACEAGORA_PERF_CONTROL_STRESS_REPEATS_FULL"
    )
    control_stress_warmup_full = _parse_nonnegative_int(
        get(ENV, "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL", "1"),
        "SPACEAGORA_PERF_CONTROL_STRESS_WARMUP_FULL"
    )
    solver_axis = _parse_solver_axis(get(ENV, "SPACEAGORA_SMART_LADDER_SOLVER_AXIS", "frozen"))
    solver_mode = lowercase(strip(get(ENV, "SPACEAGORA_SMART_LADDER_SOLVER_MODE", "auto_stiff")))
    solver_factor_modes = _parse_solver_factor_modes(get(ENV, "SPACEAGORA_SMART_LADDER_SOLVER_FACTORS", "auto_stiff,split_imex,multirate"))
    trajectory_output = _parse_bool_token(get(ENV, "SPACEAGORA_SMART_LADDER_TRAJECTORY_OUTPUT", "0"))

    for arg in ARGS
        if arg == "smoke"
            profile_name = "quick"
            ENV["SPACEAGORA_SMART_LADDER_SMOKE"] = "1"
            ENV["SPACEAGORA_PERF_SMOKE"] = "1"
            clean = true
            passes = 1
            include_layer_attribution = false
            rung_filter = String[]
            include_control_stress_per_orbit = false
            control_stress_repeats_full = 1
            control_stress_warmup_full = 0
        elseif arg in ("quick", "full")
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
        elseif startswith(arg, "--layer-attribution=")
            include_layer_attribution = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--rungs=")
            rung_filter = _parse_rung_filter(String(split(arg, "=", limit=2)[2]))
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
        elseif startswith(arg, "--solver-axis=")
            solver_axis = _parse_solver_axis(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--solver-mode=")
            solver_mode = lowercase(strip(String(split(arg, "=", limit=2)[2])))
        elseif startswith(arg, "--solver-factors=")
            solver_factor_modes = _parse_solver_factor_modes(String(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--trajectory-output=")
            trajectory_output = _parse_bool_token(String(split(arg, "=", limit=2)[2]))
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: [quick|full|smoke], --profile=..., --outdir=..., --clean=0|1, " *
                "--passes=N, --randomize-rung-order=0|1, --seed=N, --outer-only-backend=threads|process|auto, " *
                "--process-workers=N, --layer-attribution=0|1, --rungs=<csv>, --include-control-stress-per-orbit=0|1, --control-stress-repeats-full=N, " *
                "--control-stress-warmup-full=N, --solver-axis=inherit|frozen|factorial, --solver-mode=<mode>, --solver-factors=<csv>, " *
                "--trajectory-output=0|1."
            ))
        end
    end

    if isempty(solver_mode)
        throw(ArgumentError("solver-mode must be non-empty."))
    end
    if solver_axis == :factorial && isempty(solver_factor_modes)
        solver_factor_modes = [solver_mode]
    end
    if !(solver_mode in solver_factor_modes)
        pushfirst!(solver_factor_modes, solver_mode)
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
        include_layer_attribution=include_layer_attribution,
        rung_filter=rung_filter,
        include_control_stress_per_orbit=include_control_stress_per_orbit,
        control_stress_repeats_full=control_stress_repeats_full,
        control_stress_warmup_full=control_stress_warmup_full,
        solver_axis=solver_axis,
        solver_mode=solver_mode,
        solver_factor_modes=solver_factor_modes,
        trajectory_output=trajectory_output
    )
end

@inline function _ladder_matrix_modes(matrix::Symbol)::NamedTuple
    if matrix == :outer_pinned
        return (density="off", thermal="off", control="off", multibody="off", effector="off")
    elseif matrix == :full_auto
        return (density="auto", thermal="auto", control="auto", multibody="auto", effector="auto")
    elseif matrix == :attribution_density
        return (density="auto", thermal="off", control="off", multibody="off", effector="off")
    elseif matrix == :attribution_thermal
        return (density="off", thermal="auto", control="off", multibody="off", effector="off")
    elseif matrix == :attribution_control
        return (density="off", thermal="off", control="auto", multibody="off", effector="off")
    elseif matrix == :attribution_multibody
        return (density="off", thermal="off", control="off", multibody="auto", effector="off")
    elseif matrix == :attribution_effector
        return (density="off", thermal="off", control="off", multibody="off", effector="auto")
    elseif matrix == :attribution_density_thermal
        return (density="auto", thermal="auto", control="off", multibody="off", effector="off")
    elseif matrix == :attribution_density_multibody
        return (density="auto", thermal="off", control="off", multibody="auto", effector="off")
    elseif matrix == :attribution_control_effector
        return (density="off", thermal="off", control="auto", multibody="off", effector="auto")
    elseif matrix == :attribution_multibody_effector
        return (density="off", thermal="off", control="off", multibody="auto", effector="auto")
    end
    throw(ArgumentError(
        "Unsupported ladder matrix '$matrix'. Use: outer_pinned, full_auto, attribution_*."
    ))
end

@inline function _layer_attribution_backend(config::SmartLadderConfig)::String
    # Layer attribution should keep outer routing fixed so per-layer deltas are interpretable.
    return config.outer_only_backend == "auto" ? "threads" : config.outer_only_backend
end
