module Spec

using TOML

export ParameterKind, ParameterSpec, BudgetSpec, CalibrationSpec
export load_spec, default_spec, validate_spec, spec_to_dict, primary_manifest_path

@enum ParameterKind begin
    continuous
    integer
    categorical
end

Base.@kwdef struct ParameterSpec
    name::String
    kind::ParameterKind = continuous
    lower::Float64 = 0.0
    upper::Float64 = 1.0
    choices::Vector{String} = String[]
    manifest_targets::Vector{String} = String[]
    env_targets::Vector{String} = String[]
    transform::String = "set"
end

Base.@kwdef struct BudgetSpec
    initial_samples::Int = 16
    global_iters::Int = 32
    batch_size::Int = 1
    global_acquisition::String = "ei"
    bo_pool_size::Int = 256
    bo_length_scale::Float64 = 0.35
    bo_noise::Float64 = 1.0e-6
    bo_kappa::Float64 = 1.96
    bo_xi::Float64 = 0.01
    local_refine_strategy::String = "trust_region"
    local_refine_topk::Int = 4
    local_refine_steps::Int = 8
    local_refine_neighbors::Int = 4
    local_refine_init_scale::Float64 = 0.15
    local_refine_shrink::Float64 = 0.6
    local_refine_expand::Float64 = 1.2
    local_refine_min_improvement::Float64 = 1.0e-8
    robustness_samples::Int = 24
    initial_design::String = "lhs"
    robustness_uncertainty::String = "none"
    uncertainty_atm_scale::Float64 = 0.03
    uncertainty_ic_scale::Float64 = 0.002
    robust_ranking::String = "weighted"
    robust_cvar_alpha::Float64 = 0.9
    robust_p95_weight::Float64 = 0.5
    robust_fail_weight::Float64 = 5.0
    objective_huber_delta::Float64 = 1.0
    objective_lambda_fail::Float64 = 25.0
    objective_lambda_time::Float64 = 2.0
    objective_runtime_budget_s::Float64 = 1800.0
    objective_telemetry_noise_frac::Float64 = 0.01
end

Base.@kwdef struct CalibrationSpec
    schema_version::Int = 1
    id::String = "default"
    name::String = "SpaceAGORA Calibration"
    description::String = "Calibration campaign"
    output_root::String = "output/calibration"
    seed::Int = 42
    objective::String = "score"
    verification_script::String = "scripts/verify_telemetry.jl"
    manifest_paths::Vector{String} = String["test/telemetry_benchmark_manifest.toml"]
    scenario_weights::Dict{String, Float64} = Dict{String, Float64}()
    parameters::Vector{ParameterSpec} = ParameterSpec[]
    budgets::BudgetSpec = BudgetSpec()
end

@inline function primary_manifest_path(spec::CalibrationSpec)::String
    isempty(spec.manifest_paths) && throw(ArgumentError("spec.manifest_paths cannot be empty."))
    return spec.manifest_paths[1]
end

@inline function _parse_kind(raw::AbstractString)::ParameterKind
    token = lowercase(strip(String(raw)))
    if token == "continuous"
        return continuous
    elseif token == "integer"
        return integer
    elseif token == "categorical"
        return categorical
    end
    throw(ArgumentError("Unsupported parameter kind '$raw'."))
end

function _parse_parameter(tbl)::ParameterSpec
    name = String(get(tbl, "name", ""))
    isempty(name) && throw(ArgumentError("Parameter is missing non-empty field 'name'."))
    kind = _parse_kind(String(get(tbl, "kind", "continuous")))
    lower = Float64(get(tbl, "lower", 0.0))
    upper = Float64(get(tbl, "upper", 1.0))
    choices = [String(v) for v in get(tbl, "choices", Any[])]
    manifest_targets = [String(v) for v in get(tbl, "manifest_targets", Any[])]
    env_targets = [String(v) for v in get(tbl, "env_targets", Any[])]
    transform = lowercase(String(get(tbl, "transform", "set")))
    return ParameterSpec(
        name=name,
        kind=kind,
        lower=lower,
        upper=upper,
        choices=choices,
        manifest_targets=manifest_targets,
        env_targets=env_targets,
        transform=transform
    )
end

function _parse_manifest_paths(doc)::Vector{String}
    if haskey(doc, "manifest_paths")
        raw = doc["manifest_paths"]
        raw isa AbstractVector || throw(ArgumentError("manifest_paths must be an array of strings."))
        out = [String(v) for v in raw]
        return out
    elseif haskey(doc, "manifest_path")
        return [String(doc["manifest_path"])]
    end
    return String["test/telemetry_benchmark_manifest.toml"]
end

function load_spec(path::AbstractString)::CalibrationSpec
    doc = TOML.parsefile(path)

    params_raw = get(doc, "parameters", Any[])
    params = ParameterSpec[]
    for p in params_raw
        p isa AbstractDict || throw(ArgumentError("Every item in 'parameters' must be a table."))
        push!(params, _parse_parameter(p))
    end

    budgets_raw = get(doc, "budgets", Dict{String, Any}())
    budgets = BudgetSpec(
        initial_samples=Int(get(budgets_raw, "initial_samples", 16)),
        global_iters=Int(get(budgets_raw, "global_iters", 32)),
        batch_size=Int(get(budgets_raw, "batch_size", 1)),
        global_acquisition=lowercase(String(get(budgets_raw, "global_acquisition", "ei"))),
        bo_pool_size=Int(get(budgets_raw, "bo_pool_size", 256)),
        bo_length_scale=Float64(get(budgets_raw, "bo_length_scale", 0.35)),
        bo_noise=Float64(get(budgets_raw, "bo_noise", 1.0e-6)),
        bo_kappa=Float64(get(budgets_raw, "bo_kappa", 1.96)),
        bo_xi=Float64(get(budgets_raw, "bo_xi", 0.01)),
        local_refine_strategy=lowercase(String(get(budgets_raw, "local_refine_strategy", "trust_region"))),
        local_refine_topk=Int(get(budgets_raw, "local_refine_topk", 4)),
        local_refine_steps=Int(get(budgets_raw, "local_refine_steps", 8)),
        local_refine_neighbors=Int(get(budgets_raw, "local_refine_neighbors", 4)),
        local_refine_init_scale=Float64(get(budgets_raw, "local_refine_init_scale", 0.15)),
        local_refine_shrink=Float64(get(budgets_raw, "local_refine_shrink", 0.6)),
        local_refine_expand=Float64(get(budgets_raw, "local_refine_expand", 1.2)),
        local_refine_min_improvement=Float64(get(budgets_raw, "local_refine_min_improvement", 1.0e-8)),
        robustness_samples=Int(get(budgets_raw, "robustness_samples", 24)),
        initial_design=lowercase(String(get(budgets_raw, "initial_design", "lhs"))),
        robustness_uncertainty=lowercase(String(get(budgets_raw, "robustness_uncertainty", "none"))),
        uncertainty_atm_scale=Float64(get(budgets_raw, "uncertainty_atm_scale", get(budgets_raw, "robustness_jitter_scale", 0.03))),
        uncertainty_ic_scale=Float64(get(budgets_raw, "uncertainty_ic_scale", 0.002)),
        robust_ranking=lowercase(String(get(budgets_raw, "robust_ranking", "weighted"))),
        robust_cvar_alpha=Float64(get(budgets_raw, "robust_cvar_alpha", 0.9)),
        robust_p95_weight=Float64(get(budgets_raw, "robust_p95_weight", 0.5)),
        robust_fail_weight=Float64(get(budgets_raw, "robust_fail_weight", 5.0)),
        objective_huber_delta=Float64(get(budgets_raw, "objective_huber_delta", get(budgets_raw, "huber_delta", 1.0))),
        objective_lambda_fail=Float64(get(budgets_raw, "objective_lambda_fail", get(budgets_raw, "lambda_fail", 25.0))),
        objective_lambda_time=Float64(get(budgets_raw, "objective_lambda_time", get(budgets_raw, "lambda_time", 2.0))),
        objective_runtime_budget_s=Float64(get(budgets_raw, "objective_runtime_budget_s", get(budgets_raw, "runtime_budget_s", 1800.0))),
        objective_telemetry_noise_frac=Float64(get(budgets_raw, "objective_telemetry_noise_frac", get(budgets_raw, "telemetry_noise_frac", 0.01)))
    )

    weights_raw = get(doc, "scenario_weights", Dict{String, Any}())
    weights = Dict{String, Float64}(String(k) => Float64(v) for (k, v) in pairs(weights_raw))

    spec = CalibrationSpec(
        schema_version=Int(get(doc, "schema_version", get(doc, "version", 1))),
        id=String(get(doc, "id", "default")),
        name=String(get(doc, "name", "SpaceAGORA Calibration")),
        description=String(get(doc, "description", "Calibration campaign")),
        output_root=String(get(doc, "output_root", "output/calibration")),
        seed=Int(get(doc, "seed", 42)),
        objective=String(get(doc, "objective", "score")),
        verification_script=String(get(doc, "verification_script", "scripts/verify_telemetry.jl")),
        manifest_paths=_parse_manifest_paths(doc),
        scenario_weights=weights,
        parameters=params,
        budgets=budgets
    )
    validate_spec(spec)
    return spec
end

function default_spec()::CalibrationSpec
    return CalibrationSpec(
        schema_version=1,
        id="demo",
        name="Demo Calibration",
        description="Default calibration specification scaffold",
        parameters=[
            ParameterSpec(
                name="odyssey_ra_scale",
                kind=continuous,
                lower=0.995,
                upper=1.005,
                transform="mul",
                manifest_targets=["scenarios[name=odyssey].ra_m"]
            ),
            ParameterSpec(
                name="odyssey_rp_offset_m",
                kind=continuous,
                lower=-3_000.0,
                upper=3_000.0,
                transform="add",
                manifest_targets=["scenarios[name=odyssey].rp_altitude_m"]
            ),
            ParameterSpec(
                name="odyssey_gravity_model",
                kind=categorical,
                choices=["inverse_squared", "inverse_squared_j2"],
                transform="set",
                manifest_targets=["scenarios[name=odyssey].gravity_model"]
            )
        ],
        scenario_weights=Dict("odyssey" => 1.0, "vex" => 1.0, "earth_gmat" => 1.0)
    )
end

function validate_spec(spec::CalibrationSpec)::Nothing
    spec.schema_version == 1 || throw(ArgumentError(
        "Unsupported schema_version=$(spec.schema_version); expected 1."
    ))
    isempty(strip(spec.id)) && throw(ArgumentError("spec.id cannot be empty."))
    isempty(spec.parameters) && throw(ArgumentError("spec.parameters cannot be empty."))
    isempty(spec.manifest_paths) && throw(ArgumentError("spec.manifest_paths cannot be empty."))

    for p in spec.manifest_paths
        isempty(strip(p)) && throw(ArgumentError("spec.manifest_paths cannot contain empty entries."))
    end

    names_seen = Set{String}()
    for p in spec.parameters
        isempty(strip(p.name)) && throw(ArgumentError("parameter name cannot be empty."))
        p.name in names_seen && throw(ArgumentError("duplicate parameter name: $(p.name)"))
        push!(names_seen, p.name)

        if p.kind in (continuous, integer)
            p.upper >= p.lower || throw(ArgumentError("parameter $(p.name) has upper < lower."))
        elseif p.kind == categorical
            isempty(p.choices) && throw(ArgumentError("categorical parameter $(p.name) requires non-empty choices."))
        end

        p.transform in ("set", "add", "mul") || throw(ArgumentError(
            "parameter $(p.name) has unsupported transform='$(p.transform)'; use set|add|mul."
        ))
        isempty(p.manifest_targets) && isempty(p.env_targets) && throw(ArgumentError(
            "parameter $(p.name) must define at least one manifest target or env target."
        ))
    end

    b = spec.budgets
    b.initial_samples > 0 || throw(ArgumentError("budgets.initial_samples must be > 0."))
    b.global_iters >= 0 || throw(ArgumentError("budgets.global_iters must be >= 0."))
    b.batch_size > 0 || throw(ArgumentError("budgets.batch_size must be > 0."))
    b.global_acquisition in ("ei", "lcb") || throw(ArgumentError(
        "budgets.global_acquisition must be ei|lcb, got '$(b.global_acquisition)'."
    ))
    b.bo_pool_size > 0 || throw(ArgumentError("budgets.bo_pool_size must be > 0."))
    b.bo_length_scale > 0.0 || throw(ArgumentError("budgets.bo_length_scale must be > 0."))
    b.bo_noise >= 0.0 || throw(ArgumentError("budgets.bo_noise must be >= 0."))
    b.bo_kappa >= 0.0 || throw(ArgumentError("budgets.bo_kappa must be >= 0."))
    b.bo_xi >= 0.0 || throw(ArgumentError("budgets.bo_xi must be >= 0."))
    b.local_refine_strategy in ("trust_region", "perturb") || throw(ArgumentError(
        "budgets.local_refine_strategy must be trust_region|perturb, got '$(b.local_refine_strategy)'."
    ))
    b.local_refine_topk > 0 || throw(ArgumentError("budgets.local_refine_topk must be > 0."))
    b.local_refine_steps >= 0 || throw(ArgumentError("budgets.local_refine_steps must be >= 0."))
    b.local_refine_neighbors > 0 || throw(ArgumentError("budgets.local_refine_neighbors must be > 0."))
    b.local_refine_init_scale > 0.0 || throw(ArgumentError("budgets.local_refine_init_scale must be > 0."))
    b.local_refine_shrink > 0.0 || throw(ArgumentError("budgets.local_refine_shrink must be > 0."))
    b.local_refine_expand > 0.0 || throw(ArgumentError("budgets.local_refine_expand must be > 0."))
    b.local_refine_min_improvement >= 0.0 || throw(ArgumentError("budgets.local_refine_min_improvement must be >= 0."))
    b.robustness_samples >= 0 || throw(ArgumentError("budgets.robustness_samples must be >= 0."))
    b.initial_design in ("random", "lhs") || throw(ArgumentError(
        "budgets.initial_design must be random|lhs, got '$(b.initial_design)'."
    ))
    b.robustness_uncertainty in ("none", "normal") || throw(ArgumentError(
        "budgets.robustness_uncertainty must be none|normal, got '$(b.robustness_uncertainty)'."
    ))
    b.uncertainty_atm_scale >= 0.0 || throw(ArgumentError("budgets.uncertainty_atm_scale must be >= 0."))
    b.uncertainty_ic_scale >= 0.0 || throw(ArgumentError("budgets.uncertainty_ic_scale must be >= 0."))
    b.robust_ranking in ("weighted", "cvar", "p95") || throw(ArgumentError(
        "budgets.robust_ranking must be weighted|cvar|p95, got '$(b.robust_ranking)'."
    ))
    (0.0 < b.robust_cvar_alpha < 1.0) || throw(ArgumentError("budgets.robust_cvar_alpha must be in (0,1)."))
    b.robust_p95_weight >= 0.0 || throw(ArgumentError("budgets.robust_p95_weight must be >= 0."))
    b.robust_fail_weight >= 0.0 || throw(ArgumentError("budgets.robust_fail_weight must be >= 0."))
    b.objective_huber_delta > 0.0 || throw(ArgumentError("budgets.objective_huber_delta must be > 0."))
    b.objective_lambda_fail >= 0.0 || throw(ArgumentError("budgets.objective_lambda_fail must be >= 0."))
    b.objective_lambda_time >= 0.0 || throw(ArgumentError("budgets.objective_lambda_time must be >= 0."))
    b.objective_runtime_budget_s > 0.0 || throw(ArgumentError("budgets.objective_runtime_budget_s must be > 0."))
    b.objective_telemetry_noise_frac >= 0.0 || throw(ArgumentError("budgets.objective_telemetry_noise_frac must be >= 0."))
    return nothing
end

function spec_to_dict(spec::CalibrationSpec)::Dict{String, Any}
    params = Vector{Dict{String, Any}}(undef, length(spec.parameters))
    for i in eachindex(spec.parameters)
        p = spec.parameters[i]
        params[i] = Dict(
            "name" => p.name,
            "kind" => String(Symbol(p.kind)),
            "lower" => p.lower,
            "upper" => p.upper,
            "choices" => copy(p.choices),
            "manifest_targets" => copy(p.manifest_targets),
            "env_targets" => copy(p.env_targets),
            "transform" => p.transform
        )
    end

    return Dict(
        "schema_version" => spec.schema_version,
        "id" => spec.id,
        "name" => spec.name,
        "description" => spec.description,
        "output_root" => spec.output_root,
        "seed" => spec.seed,
        "objective" => spec.objective,
        "verification_script" => spec.verification_script,
        "manifest_paths" => copy(spec.manifest_paths),
        "scenario_weights" => Dict(spec.scenario_weights),
        "parameters" => params,
        "budgets" => Dict(
            "initial_samples" => spec.budgets.initial_samples,
            "global_iters" => spec.budgets.global_iters,
            "batch_size" => spec.budgets.batch_size,
            "global_acquisition" => spec.budgets.global_acquisition,
            "bo_pool_size" => spec.budgets.bo_pool_size,
            "bo_length_scale" => spec.budgets.bo_length_scale,
            "bo_noise" => spec.budgets.bo_noise,
            "bo_kappa" => spec.budgets.bo_kappa,
            "bo_xi" => spec.budgets.bo_xi,
            "local_refine_strategy" => spec.budgets.local_refine_strategy,
            "local_refine_topk" => spec.budgets.local_refine_topk,
            "local_refine_steps" => spec.budgets.local_refine_steps,
            "local_refine_neighbors" => spec.budgets.local_refine_neighbors,
            "local_refine_init_scale" => spec.budgets.local_refine_init_scale,
            "local_refine_shrink" => spec.budgets.local_refine_shrink,
            "local_refine_expand" => spec.budgets.local_refine_expand,
            "local_refine_min_improvement" => spec.budgets.local_refine_min_improvement,
            "robustness_samples" => spec.budgets.robustness_samples,
            "initial_design" => spec.budgets.initial_design,
            "robustness_uncertainty" => spec.budgets.robustness_uncertainty,
            "uncertainty_atm_scale" => spec.budgets.uncertainty_atm_scale,
            "uncertainty_ic_scale" => spec.budgets.uncertainty_ic_scale,
            "robust_ranking" => spec.budgets.robust_ranking,
            "robust_cvar_alpha" => spec.budgets.robust_cvar_alpha,
            "robust_p95_weight" => spec.budgets.robust_p95_weight,
            "robust_fail_weight" => spec.budgets.robust_fail_weight,
            "objective_huber_delta" => spec.budgets.objective_huber_delta,
            "objective_lambda_fail" => spec.budgets.objective_lambda_fail,
            "objective_lambda_time" => spec.budgets.objective_lambda_time,
            "objective_runtime_budget_s" => spec.budgets.objective_runtime_budget_s,
            "objective_telemetry_noise_frac" => spec.budgets.objective_telemetry_noise_frac
        )
    )
end

end # module Spec
