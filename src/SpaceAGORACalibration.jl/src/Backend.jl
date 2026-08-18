module Backend

using Dates
using Random
using Statistics
using TOML
using CSV
import SpaceAGORA

using ..Spec: CalibrationSpec, ParameterSpec, primary_manifest_path, continuous, integer, categorical
using ..ParamSpace: Candidate

export AbstractBackend, CommandBackend, InProcessBackend, MockBackend, BackendEvaluation
export CandidateRuntimePolicy, backend_parallel_profile, backend_full_auto_requested
export evaluate_candidate, apply_candidate_to_manifest!

abstract type AbstractBackend end

Base.@kwdef struct BackendEvaluation
    candidate_id::Int
    stage::String
    success::Bool
    score::Float64
    metrics::Dict{String, Float64} = Dict{String, Float64}()
    runtime_s::Float64 = 0.0
    error_message::String = ""
    artifacts::Dict{String, String} = Dict{String, String}()
end

Base.@kwdef struct CandidateRuntimePolicy
    outer_parallel_active::Bool = false
    outer_backend::Symbol = :none
    inner_thread_budget::Int = 1
end

Base.@kwdef struct MockBackend <: AbstractBackend
    noise_sigma::Float64 = 0.025
    fail_rate::Float64 = 0.0
end

Base.@kwdef struct CommandBackend <: AbstractBackend
    julia_cmd::Cmd = Base.julia_cmd()
    project_path::String = ".AGORA"
    verification_script::String = "src/analysis/verification/TelemetryVerification.jl"
    manifest_path::String = ""
    profile::String = "quick"
    parallel_profile::Union{Nothing, String} = nothing
    enforce::Bool = false
    plots::Union{Nothing, Bool} = nothing
end

Base.@kwdef struct InProcessBackend{F} <: AbstractBackend
    run_verification::F = SpaceAGORA.run_verification
    manifest_path::String = ""
    profile::String = "quick"
    parallel_profile::Union{Nothing, String} = nothing
    enforce::Bool = false
    plots::Union{Nothing, Bool} = nothing
end

@inline function _profile_for_stage(stage::String, backend::CommandBackend)::String
    if stage in ("local_refine_full", "robustness_validation", "promote")
        return "full"
    end
    return backend.profile
end

@inline function _profile_for_stage(stage::String, backend::InProcessBackend)::String
    if stage in ("local_refine_full", "robustness_validation", "promote")
        return "full"
    end
    return backend.profile
end

@inline function _replica_id(candidate::Candidate)::Union{Nothing, Int}
    if !haskey(candidate.values, "robustness_replica")
        return nothing
    end
    raw = candidate.values["robustness_replica"]
    return try
        Int(raw)
    catch
        nothing
    end
end

@inline function _lookup_parameter(spec::CalibrationSpec, name::String)::ParameterSpec
    for p in spec.parameters
        if p.name == name
            return p
        end
    end
    throw(ArgumentError("Candidate parameter '$name' not found in spec.parameters."))
end

@inline function _coerce_parameter_value(param::ParameterSpec, raw)
    if param.kind == continuous
        return Float64(raw)
    elseif param.kind == integer
        return Int(raw)
    elseif param.kind == categorical
        val = String(raw)
        val in param.choices || throw(ArgumentError(
            "Categorical parameter $(param.name) value '$val' is not in choices=$(param.choices)."
        ))
        return val
    end
    throw(ArgumentError("Unsupported parameter kind for $(param.name)."))
end

@inline function _parse_segment(seg::AbstractString)
    s = strip(seg)
    isempty(s) && throw(ArgumentError("Manifest target contains an empty path segment."))

    m = match(r"^([A-Za-z0-9_-]+)(?:\[(.+)\])?$", s)
    m === nothing && throw(ArgumentError("Invalid target segment '$seg'."))
    key = String(m.captures[1])
    selector_raw = m.captures[2]
    if selector_raw === nothing
        return key, nothing
    end

    selector = strip(String(selector_raw))
    if selector == "*"
        return key, (:all, nothing)
    end
    if occursin("=", selector)
        parts = split(selector, "=", limit=2)
        field = strip(parts[1])
        value = strip(parts[2])
        isempty(field) && throw(ArgumentError("Invalid selector field in '$seg'."))
        return key, (:field_eq, (field, value))
    end

    idx = try
        parse(Int, selector)
    catch
        throw(ArgumentError("Invalid selector '$selector' in target segment '$seg'."))
    end
    idx >= 1 || throw(ArgumentError("Selector index in '$seg' must be >= 1."))
    return key, (:index, idx)
end

function _select_items(items, selector, seg::AbstractString)
    mode, payload = selector
    if mode == :all
        return collect(items)
    elseif mode == :index
        idx = Int(payload)
        idx <= length(items) || throw(ArgumentError(
            "Selector index out of bounds in '$seg': index=$idx length=$(length(items))."
        ))
        return Any[items[idx]]
    elseif mode == :field_eq
        field, expected = payload
        out = Any[]
        for item in items
            item isa AbstractDict || continue
            if haskey(item, field) && string(item[field]) == expected
                push!(out, item)
            end
        end
        isempty(out) && throw(ArgumentError(
            "Selector '$field=$expected' in '$seg' matched zero entries."
        ))
        return out
    end
    throw(ArgumentError("Unsupported selector mode '$mode'."))
end

function _resolve_target_parents(doc::AbstractDict, target::String)
    parts = split(target, '.')
    length(parts) >= 1 || throw(ArgumentError("Target path '$target' is invalid."))

    nodes = Any[doc]
    for seg in parts[1:end-1]
        key, selector = _parse_segment(seg)
        next_nodes = Any[]
        for node in nodes
            node isa AbstractDict || throw(ArgumentError("Path segment '$seg' expected a table parent."))
            haskey(node, key) || throw(ArgumentError("Path segment '$seg' key '$key' missing in target '$target'."))
            child = node[key]
            if selector === nothing
                push!(next_nodes, child)
            else
                child isa AbstractVector || throw(ArgumentError(
                    "Path segment '$seg' selector requires an array in target '$target'."
                ))
                append!(next_nodes, _select_items(child, selector, seg))
            end
        end
        nodes = next_nodes
    end

    final_key, final_selector = _parse_segment(parts[end])
    final_selector === nothing || throw(ArgumentError(
        "Final target segment in '$target' must be a plain key (no selector)."
    ))

    out = Tuple{AbstractDict, String}[]
    for node in nodes
        node isa AbstractDict || throw(ArgumentError("Target '$target' resolved to non-table parent."))
        push!(out, (node, final_key))
    end
    return out
end

@inline function _apply_transform(existing, value, transform::String)
    if transform == "set"
        return value
    elseif transform == "add"
        if existing isa AbstractVector
            return [Float64(x) + Float64(value) for x in existing]
        end
        return Float64(existing) + Float64(value)
    elseif transform == "mul"
        if existing isa AbstractVector
            return [Float64(x) * Float64(value) for x in existing]
        end
        return Float64(existing) * Float64(value)
    end
    throw(ArgumentError("Unsupported transform '$transform'."))
end

function _apply_candidate_to_manifest!(
    manifest_doc::Dict{String, Any},
    spec::CalibrationSpec,
    candidate::Candidate
)::Nothing
    for (pname, raw_value) in candidate.values
        pname == "robustness_replica" && continue
        param = _lookup_parameter(spec, pname)
        value = _coerce_parameter_value(param, raw_value)

        for target in param.manifest_targets
            parents = _resolve_target_parents(manifest_doc, target)
            for (parent, key) in parents
                if param.transform == "set"
                    parent[key] = value
                else
                    haskey(parent, key) || throw(ArgumentError(
                        "Target '$target' requires existing numeric value for transform=$(param.transform)."
                    ))
                    parent[key] = _apply_transform(parent[key], value, param.transform)
                end
            end
        end
    end
    return nothing
end

function apply_candidate_to_manifest!(
    manifest_doc::Dict{String, Any},
    spec::CalibrationSpec,
    candidate::Candidate
)::Nothing
    return _apply_candidate_to_manifest!(manifest_doc, spec, candidate)
end

@inline function _f64(v)::Float64
    if v isa Integer
        return Float64(v)
    elseif v isa AbstractFloat
        return Float64(v)
    end
    return Float64(hash(string(v)) % 1000) / 1000.0
end

function _mock_score(c::Candidate)::Float64
    vals = collect(values(c.values))
    isempty(vals) && return Inf
    fvals = [_f64(v) for v in vals]
    center = mean(fvals)
    spread = std(fvals)
    return abs(center - 0.42) + 0.25 * spread
end

function evaluate_candidate(
    backend::MockBackend,
    spec::CalibrationSpec,
    candidate::Candidate;
    stage::String=candidate.stage,
    run_dir::Union{Nothing, String}=nothing,
    runtime_policy::Union{Nothing, CandidateRuntimePolicy}=nothing
)::BackendEvaluation
    _ = run_dir
    _ = runtime_policy
    rng = MersenneTwister(hash((spec.seed, candidate.id, stage)))
    if rand(rng) < backend.fail_rate
        return BackendEvaluation(
            candidate_id=candidate.id,
            stage=stage,
            success=false,
            score=Inf,
            runtime_s=0.0,
            error_message="mock_failure"
        )
    end

    noise = backend.noise_sigma > 0 ? backend.noise_sigma * randn(rng) : 0.0
    score = max(_mock_score(candidate) + noise, 0.0)
    return BackendEvaluation(
        candidate_id=candidate.id,
        stage=stage,
        success=true,
        score=score,
        runtime_s=0.0,
        metrics=Dict("mock_score" => score)
    )
end

@inline function _safe_float(v; default::Float64=NaN)::Float64
    if v === missing || v === nothing
        return default
    elseif v isa AbstractFloat || v isa Integer
        return Float64(v)
    elseif v isa AbstractString
        t = strip(String(v))
        isempty(t) && return default
        parsed = try
            parse(Float64, t)
        catch
            default
        end
        return parsed
    end
    return default
end

@inline function _safe_bool(v; default::Bool=false)::Bool
    if v === missing || v === nothing
        return default
    elseif v isa Bool
        return v
    elseif v isa Integer
        return v != 0
    elseif v isa AbstractString
        token = lowercase(strip(String(v)))
        if token in ("1", "true", "yes", "on")
            return true
        elseif token in ("0", "false", "no", "off")
            return false
        end
    end
    return default
end

@inline function _huber(x::Float64, delta::Float64)::Float64
    ax = abs(x)
    if ax <= delta
        return 0.5 * x * x
    end
    return delta * (ax - 0.5 * delta)
end

function _runtime_from_rows(rows, wall_runtime_s::Float64)::Float64
    for row in rows
        if hasproperty(row, :total_runtime_s)
            v = _safe_float(getproperty(row, :total_runtime_s))
            if isfinite(v) && v > 0.0
                return v
            end
        end
    end
    for row in rows
        if hasproperty(row, :simulation_runtime_s)
            v = _safe_float(getproperty(row, :simulation_runtime_s))
            if isfinite(v) && v > 0.0
                return v
            end
        end
    end
    return wall_runtime_s
end

@inline function _summary_supports_robust_objective(rows)::Bool
    isempty(rows) && return false
    row = rows[1]
    return hasproperty(row, :scenario) &&
           hasproperty(row, :rmse_km) &&
           hasproperty(row, :max_abs_km) &&
           hasproperty(row, :limit_max_abs_km) &&
           (hasproperty(row, :limit_max_rmse_km) || hasproperty(row, :limit_nmae))
end

@inline function _row_metric_value(
    row,
    legacy_name::Symbol,
    display_name::Symbol
)::Float64
    if hasproperty(row, display_name)
        value = _safe_float(getproperty(row, display_name))
        if isfinite(value)
            return value
        end
    end
    if hasproperty(row, legacy_name)
        return _safe_float(getproperty(row, legacy_name))
    end
    return NaN
end

function _objective_from_rows(
    rows,
    spec::CalibrationSpec;
    run_failed::Bool,
    runtime_s::Float64,
    noise_rng::Union{Nothing, AbstractRNG}=nothing
)::NamedTuple
    budgets = spec.budgets
    _ = runtime_s

    if run_failed || isempty(rows)
        return (
            score=Inf,
            base_loss=Inf,
            fail_penalty=0.0,
            runtime_penalty=0.0,
            all_pass=false,
            failed_rows=0.0
        )
    end

    if !_summary_supports_robust_objective(rows)
        nmae_values = Float64[]
        failed_rows = 0.0
        all_pass = true
        for row in rows
            if hasproperty(row, :nmae)
                nmae = _safe_float(getproperty(row, :nmae))
                isfinite(nmae) && push!(nmae_values, nmae)
            end
            if hasproperty(row, :pass)
                row_pass = _safe_bool(getproperty(row, :pass); default=true)
                all_pass &= row_pass
                !row_pass && (failed_rows += 1.0)
            end
        end

        base_loss = isempty(nmae_values) ? Inf : mean(nmae_values)
        return (
            score=base_loss,
            base_loss=base_loss,
            fail_penalty=0.0,
            runtime_penalty=0.0,
            all_pass=all_pass,
            failed_rows=failed_rows
        )
    end

    base_loss = 0.0
    all_pass = true
    failed_rows = 0.0

    for row in rows
        scenario = hasproperty(row, :scenario) ? String(getproperty(row, :scenario)) : "default"
        weight = get(spec.scenario_weights, scenario, 1.0)

        # Prefer display-normalized metrics when available (e.g., velocity rows in m/s).
        rmse = _row_metric_value(row, :rmse_km, :rmse_display)
        max_abs = _row_metric_value(row, :max_abs_km, :max_abs_display)
        nmae = hasproperty(row, :nmae) ? _safe_float(getproperty(row, :nmae)) : NaN
        lim_rmse = _row_metric_value(row, :limit_max_rmse_km, :limit_max_rmse_display)
        lim_abs = _row_metric_value(row, :limit_max_abs_km, :limit_max_abs_display)
        lim_nmae = hasproperty(row, :limit_nmae) ? _safe_float(getproperty(row, :limit_nmae)) : NaN

        rmse_ratio = if isfinite(lim_rmse) && lim_rmse > 0.0
            rmse / lim_rmse
        elseif isfinite(lim_nmae) && lim_nmae > 0.0 && isfinite(nmae)
            nmae / lim_nmae
        else
            rmse
        end

        abs_ratio = if isfinite(lim_abs) && lim_abs > 0.0
            max_abs / lim_abs
        else
            max_abs
        end

        if noise_rng !== nothing && budgets.objective_telemetry_noise_frac > 0.0
            frac = budgets.objective_telemetry_noise_frac
            rmse_ratio = max(0.0, rmse_ratio * (1.0 + frac * randn(noise_rng)))
            abs_ratio = max(0.0, abs_ratio * (1.0 + frac * randn(noise_rng)))
        end

        term = weight * (_huber(rmse_ratio, budgets.objective_huber_delta) + 0.5 * _huber(abs_ratio, budgets.objective_huber_delta))
        base_loss += term

        if hasproperty(row, :pass)
            row_pass = _safe_bool(getproperty(row, :pass); default=true)
            all_pass &= row_pass
            !row_pass && (failed_rows += 1.0)
        end
    end

    score = base_loss

    return (
        score=score,
        base_loss=base_loss,
        fail_penalty=0.0,
        runtime_penalty=0.0,
        all_pass=all_pass,
        failed_rows=failed_rows
    )
end

@inline _wrap_0_360(deg::Float64) = mod(mod(deg, 360.0) + 360.0, 360.0)
@inline _clamp_positive(x::Float64) = max(1.0e-9, x)

function _apply_uncertainty_to_scenario!(
    sc::Dict{String, Any},
    rng::AbstractRNG,
    spec::CalibrationSpec
)::Nothing
    atm_scale = spec.budgets.uncertainty_atm_scale
    ic_scale = spec.budgets.uncertainty_ic_scale

    if haskey(sc, "atmosphere_truth")
        at = Dict{String, Any}(sc["atmosphere_truth"])
        base_seed = haskey(at, "gram_seed") ? Int(at["gram_seed"]) : 1001
        at["gram_seed"] = base_seed + rand(rng, 1:1_000_000)

        base_scales = if haskey(at, "gram_perturbation_scales")
            [Float64(v) for v in at["gram_perturbation_scales"]]
        else
            [0.0, 0.0, 0.0, 0.0]
        end
        perturbed = Float64[]
        for s in base_scales
            sigma = s > 0.0 ? atm_scale * abs(s) : atm_scale
            push!(perturbed, max(0.0, s + sigma * randn(rng)))
        end
        at["gram_perturbation_scales"] = perturbed
        sc["atmosphere_truth"] = at
    end

    if haskey(sc, "ra_m")
        sc["ra_m"] = _clamp_positive(Float64(sc["ra_m"]) * (1.0 + ic_scale * randn(rng)))
    end
    if haskey(sc, "rp_altitude_m")
        sc["rp_altitude_m"] = _clamp_positive(Float64(sc["rp_altitude_m"]) + 1500.0 * ic_scale * randn(rng))
    end
    if haskey(sc, "i_deg")
        sc["i_deg"] = clamp(Float64(sc["i_deg"]) + 0.5 * ic_scale * randn(rng), 0.0, 180.0)
    end
    if haskey(sc, "aop_deg")
        sc["aop_deg"] = _wrap_0_360(Float64(sc["aop_deg"]) + 2.0 * ic_scale * randn(rng))
    end
    if haskey(sc, "raan_deg")
        sc["raan_deg"] = _wrap_0_360(Float64(sc["raan_deg"]) + 2.0 * ic_scale * randn(rng))
    end
    if haskey(sc, "ta_deg")
        sc["ta_deg"] = _wrap_0_360(Float64(sc["ta_deg"]) + 5.0 * ic_scale * randn(rng))
    end

    return nothing
end

function _apply_uncertainty_to_manifest!(
    manifest_doc::Dict{String, Any},
    spec::CalibrationSpec,
    candidate_id::Int,
    replica::Int
)::Nothing
    spec.budgets.robustness_uncertainty == "normal" || return nothing

    if !haskey(manifest_doc, "scenarios") || !(manifest_doc["scenarios"] isa AbstractVector)
        return nothing
    end

    scenarios = manifest_doc["scenarios"]
    for i in eachindex(scenarios)
        scenarios[i] isa AbstractDict || continue
        sc = Dict{String, Any}(scenarios[i])
        sc_name = String(get(sc, "name", "scenario"))
        rng = MersenneTwister(hash((spec.seed, "robustness_uncertainty", sc_name, candidate_id, replica)))
        _apply_uncertainty_to_scenario!(sc, rng, spec)
        scenarios[i] = sc
    end
    manifest_doc["scenarios"] = scenarios
    return nothing
end

@inline function _format_env_value(v)::String
    if v isa Bool
        return v ? "1" : "0"
    elseif v isa Integer
        return string(Int(v))
    elseif v isa AbstractFloat
        return string(Float64(v))
    end
    return string(v)
end

@inline function _env_with_default(name::String, default::String)::String
    raw = strip(get(ENV, name, ""))
    return isempty(raw) ? default : raw
end

@inline function _backend_parallel_profile_raw(::AbstractBackend)::Union{Nothing, String}
    return nothing
end

@inline function _backend_parallel_profile_raw(backend::CommandBackend)::Union{Nothing, String}
    return backend.parallel_profile
end

@inline function _backend_parallel_profile_raw(backend::InProcessBackend)::Union{Nothing, String}
    return backend.parallel_profile
end

@inline function _is_full_auto_profile_token(raw::AbstractString)::Bool
    token = lowercase(strip(String(raw)))
    token = replace(token, "-" => "_")
    token = replace(token, " " => "")
    return token in (
        "r5",
        "r5_full_auto",
        "r5fullauto",
        "r5_calibration_full_auto",
        "full_smart",
        # Legacy aliases:
        "r4_full_auto",
        "r4fullauto",
        "r4_calibration_full_auto"
    )
end

@inline function backend_full_auto_requested(backend::AbstractBackend)::Bool
    raw = _backend_parallel_profile_raw(backend)
    if raw === nothing || isempty(strip(raw))
        env_raw = strip(get(ENV, "SPACEAGORA_CALIBRATION_PARALLEL_PROFILE", ""))
        return _is_full_auto_profile_token(env_raw)
    end
    return _is_full_auto_profile_token(raw)
end

function _backend_parallel_profile(backend::AbstractBackend)::SpaceAGORA.ParallelProfile
    raw = _backend_parallel_profile_raw(backend)
    if raw === nothing || isempty(strip(raw))
        env_raw = strip(get(ENV, "SPACEAGORA_CALIBRATION_PARALLEL_PROFILE", ""))
        raw = isempty(env_raw) ? "R0" : env_raw
    end
    try
        return SpaceAGORA.parse_parallel_profile(raw)
    catch err
        if _is_full_auto_profile_token(raw)
            # Compatibility fallback across SpaceAGORA versions with different full-auto tokens.
            for fallback in ("R5", "R4_full_auto", "R4")
                try
                    return SpaceAGORA.parse_parallel_profile(fallback)
                catch
                end
            end
        end
        rethrow(err)
    end
end

@inline backend_parallel_profile(backend::AbstractBackend) = _backend_parallel_profile(backend)

@inline function _calibration_parallel_evaluations(spec::CalibrationSpec)::Int
    default = max(1, spec.budgets.parallel_evaluations)
    raw = strip(get(ENV, "SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS", ""))
    isempty(raw) && return default
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS must be a positive integer, got '$raw'"))
    end
    parsed > 0 || throw(ArgumentError(
        "SPACEAGORA_CALIBRATION_PARALLEL_EVALUATIONS must be a positive integer, got $parsed"
    ))
    return parsed
end

@inline function _coerce_outer_backend_token(backend::Symbol)::String
    if backend == :none
        return "none"
    elseif backend == :threads
        return "threads"
    elseif backend == :process
        return "process"
    elseif backend == :auto
        return "auto"
    end
    throw(ArgumentError("Unsupported runtime outer backend '$backend'."))
end

function _candidate_runtime_policy_env_pairs(
    backend::AbstractBackend,
    spec::CalibrationSpec;
    runtime_policy::Union{Nothing, CandidateRuntimePolicy}=nothing
)::Vector{Pair{String, String}}
    profile = _backend_parallel_profile(backend)
    outer_active = runtime_policy === nothing ? (_calibration_parallel_evaluations(spec) > 1) : runtime_policy.outer_parallel_active
    profile_pairs = SpaceAGORA.profile_env_pairs(
        profile;
        preserve_existing=true,
        outer_parallel_active=outer_active
    )

    env_map = Dict{String, String}()
    env_map["OPENBLAS_NUM_THREADS"] = _env_with_default("OPENBLAS_NUM_THREADS", "1")
    env_map["SPACEAGORA_INNER_THREAD_BUDGET"] = if runtime_policy === nothing
        _env_with_default("SPACEAGORA_INNER_THREAD_BUDGET", "1")
    else
        string(max(1, runtime_policy.inner_thread_budget))
    end
    for (k, v) in profile_pairs
        env_map[k] = v
    end
    if runtime_policy !== nothing
        env_map["SPACEAGORA_PERF_PARALLEL_BACKEND"] = _coerce_outer_backend_token(runtime_policy.outer_backend)
        env_map["SPACEAGORA_OUTER_PARALLEL_ACTIVE"] = runtime_policy.outer_parallel_active ? "1" : "0"
    end
    return collect(env_map)
end

function _candidate_env_pairs(
    backend::AbstractBackend,
    spec::CalibrationSpec,
    candidate::Candidate;
    runtime_policy::Union{Nothing, CandidateRuntimePolicy}=nothing
)::Vector{Pair{String, String}}
    env_map = Dict{String, String}()

    for (pname, raw_value) in candidate.values
        pname == "robustness_replica" && continue
        param = _lookup_parameter(spec, pname)
        isempty(param.env_targets) && continue
        value = _coerce_parameter_value(param, raw_value)
        for env_name in param.env_targets
            env_map[env_name] = _format_env_value(value)
        end
    end

    for (k, v) in _candidate_runtime_policy_env_pairs(backend, spec; runtime_policy=runtime_policy)
        get!(env_map, k, v)
    end

    return collect(env_map)
end

@inline function _candidate_runtime_policy_env_pairs()::Vector{Pair{String, String}}
    return _candidate_runtime_policy_env_pairs(CommandBackend(), CalibrationSpec())
end

@inline function _candidate_env_pairs(spec::CalibrationSpec, candidate::Candidate)::Vector{Pair{String, String}}
    return _candidate_env_pairs(CommandBackend(), spec, candidate)
end

@inline function _profile_symbol(profile::String)::Symbol
    token = lowercase(strip(profile))
    token in ("quick", "full") || throw(ArgumentError("backend profile must be quick|full, got '$profile'."))
    return Symbol(token)
end

function _verification_request(
    backend::InProcessBackend,
    profile::String,
    tuned_manifest_path::String,
    summary_path::String,
    errors_path::String
)
    profile_symbol = _profile_symbol(profile)
    if backend.plots === nothing
        return SpaceAGORA.VerificationRequest(
            profile=profile_symbol,
            out_summary=summary_path,
            out_errors=errors_path,
            manifest_path=tuned_manifest_path,
            enforce=backend.enforce
        )
    end
    return SpaceAGORA.VerificationRequest(
        profile=profile_symbol,
        out_summary=summary_path,
        out_errors=errors_path,
        manifest_path=tuned_manifest_path,
        enforce=backend.enforce,
        generate_plots=backend.plots
    )
end

function _read_summary_rows(path::String)
    if !isfile(path)
        return Any[]
    end
    return collect(CSV.File(path; silencewarnings=true))
end

function evaluate_candidate(
    backend::InProcessBackend,
    spec::CalibrationSpec,
    candidate::Candidate;
    stage::String=candidate.stage,
    run_dir::Union{Nothing, String}=nothing,
    runtime_policy::Union{Nothing, CandidateRuntimePolicy}=nothing
)::BackendEvaluation
    workdir = run_dir === nothing ? mktempdir() : run_dir
    replica = _replica_id(candidate)
    suffix = replica === nothing ? "" : "_r$(lpad(string(replica), 3, '0'))"
    cdir = joinpath(workdir, "candidate_$(lpad(string(candidate.id), 4, '0'))_$(stage)$(suffix)")
    mkpath(cdir)

    payload_path = joinpath(cdir, "candidate.toml")
    open(payload_path, "w") do io
        TOML.print(io, Dict(
            "candidate_id" => candidate.id,
            "stage" => stage,
            "values" => Dict(candidate.values)
        ))
    end

    summary_path = joinpath(cdir, "summary.csv")
    errors_path = joinpath(cdir, "errors.csv")

    manifest_base = isempty(strip(backend.manifest_path)) ? primary_manifest_path(spec) : backend.manifest_path
    isfile(manifest_base) || throw(ArgumentError("Base manifest file not found: $manifest_base"))

    manifest_doc = TOML.parsefile(manifest_base)
    tuned_manifest_path = joinpath(cdir, "manifest_tuned.toml")
    _apply_candidate_to_manifest!(manifest_doc, spec, candidate)

    if replica !== nothing
        _apply_uncertainty_to_manifest!(manifest_doc, spec, candidate.id, replica)
    end

    open(tuned_manifest_path, "w") do io
        TOML.print(io, manifest_doc)
    end

    profile = _profile_for_stage(stage, backend)
    request = _verification_request(
        backend,
        profile,
        tuned_manifest_path,
        summary_path,
        errors_path
    )
    env_pairs = _candidate_env_pairs(backend, spec, candidate; runtime_policy=runtime_policy)

    ok = true
    err = ""
    wall_runtime_s = @elapsed begin
        try
            withenv(env_pairs...) do
                backend.run_verification(request)
            end
        catch e
            ok = false
            err = sprint(showerror, e)
        end
    end

    rows = ok ? _read_summary_rows(summary_path) : Any[]
    runtime_s = _runtime_from_rows(rows, wall_runtime_s)

    obj_rng = if replica === nothing
        nothing
    else
        MersenneTwister(hash((spec.seed, "telemetry_noise", candidate.id, replica, stage)))
    end

    obj = _objective_from_rows(
        rows,
        spec;
        run_failed=(!ok || isempty(rows)),
        runtime_s=runtime_s,
        noise_rng=obj_rng
    )

    success = ok && !isempty(rows) && isfinite(obj.score)

    return BackendEvaluation(
        candidate_id=candidate.id,
        stage=stage,
        success=success,
        score=obj.score,
        runtime_s=runtime_s,
        metrics=Dict(
            "objective_base" => obj.base_loss,
            "penalty_fail" => obj.fail_penalty,
            "penalty_time" => obj.runtime_penalty,
            "all_pass" => obj.all_pass ? 1.0 : 0.0,
            "failed_rows" => obj.failed_rows,
            "summary_rows" => Float64(length(rows))
        ),
        error_message=ok ? (isempty(rows) ? "summary_missing_or_invalid" : "") : (isempty(err) ? "inprocess_failed" : err),
        artifacts=Dict(
            "summary" => summary_path,
            "errors" => errors_path,
            "candidate" => payload_path,
            "manifest_base" => manifest_base,
            "manifest_tuned" => tuned_manifest_path
        )
    )
end

function evaluate_candidate(
    backend::CommandBackend,
    spec::CalibrationSpec,
    candidate::Candidate;
    stage::String=candidate.stage,
    run_dir::Union{Nothing, String}=nothing,
    runtime_policy::Union{Nothing, CandidateRuntimePolicy}=nothing
)::BackendEvaluation
    workdir = run_dir === nothing ? mktempdir() : run_dir
    replica = _replica_id(candidate)
    suffix = replica === nothing ? "" : "_r$(lpad(string(replica), 3, '0'))"
    cdir = joinpath(workdir, "candidate_$(lpad(string(candidate.id), 4, '0'))_$(stage)$(suffix)")
    mkpath(cdir)

    payload_path = joinpath(cdir, "candidate.toml")
    open(payload_path, "w") do io
        TOML.print(io, Dict(
            "candidate_id" => candidate.id,
            "stage" => stage,
            "values" => Dict(candidate.values)
        ))
    end

    summary_path = joinpath(cdir, "summary.csv")
    errors_path = joinpath(cdir, "errors.csv")
    log_path = joinpath(cdir, "run.log")

    manifest_base = isempty(strip(backend.manifest_path)) ? primary_manifest_path(spec) : backend.manifest_path
    isfile(manifest_base) || throw(ArgumentError("Base manifest file not found: $manifest_base"))

    manifest_doc = TOML.parsefile(manifest_base)
    tuned_manifest_path = joinpath(cdir, "manifest_tuned.toml")
    _apply_candidate_to_manifest!(manifest_doc, spec, candidate)

    if replica !== nothing
        _apply_uncertainty_to_manifest!(manifest_doc, spec, candidate.id, replica)
    end

    open(tuned_manifest_path, "w") do io
        TOML.print(io, manifest_doc)
    end

    profile = _profile_for_stage(stage, backend)
    cmd = if backend.plots === nothing
        `$(backend.julia_cmd) --project=$(backend.project_path) --startup-file=no $(backend.verification_script) --profile=$(profile) --manifest=$(tuned_manifest_path) --out-summary=$(summary_path) --out-errors=$(errors_path) --enforce=$(backend.enforce ? "1" : "0")`
    else
        `$(backend.julia_cmd) --project=$(backend.project_path) --startup-file=no $(backend.verification_script) --profile=$(profile) --manifest=$(tuned_manifest_path) --out-summary=$(summary_path) --out-errors=$(errors_path) --enforce=$(backend.enforce ? "1" : "0") --plots=$(backend.plots ? "1" : "0")`
    end

    env_pairs = _candidate_env_pairs(backend, spec, candidate; runtime_policy=runtime_policy)

    ok = true
    err = ""
    wall_runtime_s = @elapsed begin
        open(log_path, "w") do io
            try
                withenv(env_pairs...) do
                    run(pipeline(cmd, stdout=io, stderr=io))
                end
            catch e
                ok = false
                err = sprint(showerror, e)
            end
        end
    end

    rows = ok ? _read_summary_rows(summary_path) : Any[]
    runtime_s = _runtime_from_rows(rows, wall_runtime_s)

    obj_rng = if replica === nothing
        nothing
    else
        MersenneTwister(hash((spec.seed, "telemetry_noise", candidate.id, replica, stage)))
    end

    obj = _objective_from_rows(
        rows,
        spec;
        run_failed=(!ok || isempty(rows)),
        runtime_s=runtime_s,
        noise_rng=obj_rng
    )

    success = ok && !isempty(rows) && isfinite(obj.score)

    return BackendEvaluation(
        candidate_id=candidate.id,
        stage=stage,
        success=success,
        score=obj.score,
        runtime_s=runtime_s,
        metrics=Dict(
            "objective_base" => obj.base_loss,
            "penalty_fail" => obj.fail_penalty,
            "penalty_time" => obj.runtime_penalty,
            "all_pass" => obj.all_pass ? 1.0 : 0.0,
            "failed_rows" => obj.failed_rows,
            "summary_rows" => Float64(length(rows))
        ),
        error_message=ok ? (isempty(rows) ? "summary_missing_or_invalid" : "") : (isempty(err) ? "command_failed" : err),
        artifacts=Dict(
            "log" => log_path,
            "summary" => summary_path,
            "errors" => errors_path,
            "candidate" => payload_path,
            "manifest_base" => manifest_base,
            "manifest_tuned" => tuned_manifest_path
        )
    )
end

end # module Backend
