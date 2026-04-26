const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_BASE_MANIFEST = joinpath(REPO_ROOT, "test", "telemetry_benchmark_manifest.toml")
const DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "telemetry_tuning", "hybrid")
const TELEMETRY_STUDY_SCRIPT = joinpath(REPO_ROOT, "benchmarks", "studies", "telemetry_orbit_accuracy_study.jl")
const TUNER_HEARTBEAT_SECONDS = 120.0

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using Random
using Statistics
using TOML

const GRAVITY_VARIANTS = ["harm_deg2", "harm_deg4", "harm_deg8", "harm_deg20", "j2_only"]

const PARAMETER_ORDER = [
    "epoch_shift_s",
    "ra_scale",
    "rp_altitude_offset_m",
    "i_offset_deg",
    "aop_offset_deg",
    "raan_offset_deg",
    "ta_offset_deg",
    "bus_mass_scale",
    "prop_mass_scale",
    "panel_mass_scale",
    "bus_dims_scale",
    "panel_dims_scale",
    "panel_offset_scale",
    "srp_cr_scale",
    "srp_area_scale",
    "gravity_variant",
    "dv_global_scale",
    "dv_early_scale",
    "dv_orbit7_bias_mps",
    "solver_mode",
    "dt_max_orbit_s",
    "dt_max_atm_s"
]

Base.@kwdef struct TunerConfig
    base_manifest::String = DEFAULT_BASE_MANIFEST
    outdir::String = DEFAULT_OUTDIR
    n_init::Int = 16
    n_bo_iters::Int = 12
    batch_size::Int = 2
    bo_pool_size::Int = 256
    acquisition::Symbol = :lcb
    kappa_start::Float64 = 2.5
    kappa_end::Float64 = 0.8
    ei_xi::Float64 = 0.01
    surrogate_length_scale::Float64 = 0.35
    surrogate_noise::Float64 = 1.0e-6
    local_seed_topk::Int = 4
    local_steps::Int = 6
    local_step_init::Float64 = 0.12
    local_step_shrink::Float64 = 0.6
    local_step_expand::Float64 = 1.2
    local_min_improve::Float64 = 1.0e-6
    finalists::Int = 4
    robust_draws::Int = 12
    robust_alpha::Float64 = 0.5
    robust_beta::Float64 = 5.0
    uncertainty_atm_scale::Float64 = 0.03
    uncertainty_ic_scale::Float64 = 0.002
    telemetry_noise_frac::Float64 = 0.01
    huber_delta::Float64 = 1.0
    lambda_fail::Float64 = 25.0
    lambda_time::Float64 = 2.0
    runtime_budget_s::Float64 = 1800.0
    parallel_candidates::Int = 1
    seed::Int = 42
    enforce::Bool = false
    calibration_enabled::Bool = false
    init_design::Symbol = :lhs
    full_confirm::Bool = true
    solver_modes::Vector{String} = String["auto_stiff"]
    dt_max_orbit_values::Vector{Float64} = [30.0, 60.0]
    dt_max_atm_values::Vector{Float64} = [0.1, 0.2]
    maxiters_quick::Union{Nothing, Int} = nothing
    maxiters_full::Union{Nothing, Int} = nothing
    scenario_weights::Dict{String, Float64} = Dict(
        "odyssey" => 1.0,
        "vex" => 1.0,
        "earth_gmat" => 1.0
    )
end

@inline function _tuner_config_fields(cfg::TunerConfig)
    return (
        base_manifest=cfg.base_manifest,
        outdir=cfg.outdir,
        n_init=cfg.n_init,
        n_bo_iters=cfg.n_bo_iters,
        batch_size=cfg.batch_size,
        bo_pool_size=cfg.bo_pool_size,
        acquisition=cfg.acquisition,
        kappa_start=cfg.kappa_start,
        kappa_end=cfg.kappa_end,
        ei_xi=cfg.ei_xi,
        surrogate_length_scale=cfg.surrogate_length_scale,
        surrogate_noise=cfg.surrogate_noise,
        local_seed_topk=cfg.local_seed_topk,
        local_steps=cfg.local_steps,
        local_step_init=cfg.local_step_init,
        local_step_shrink=cfg.local_step_shrink,
        local_step_expand=cfg.local_step_expand,
        local_min_improve=cfg.local_min_improve,
        finalists=cfg.finalists,
        robust_draws=cfg.robust_draws,
        robust_alpha=cfg.robust_alpha,
        robust_beta=cfg.robust_beta,
        uncertainty_atm_scale=cfg.uncertainty_atm_scale,
        uncertainty_ic_scale=cfg.uncertainty_ic_scale,
        telemetry_noise_frac=cfg.telemetry_noise_frac,
        huber_delta=cfg.huber_delta,
        lambda_fail=cfg.lambda_fail,
        lambda_time=cfg.lambda_time,
        runtime_budget_s=cfg.runtime_budget_s,
        parallel_candidates=cfg.parallel_candidates,
        seed=cfg.seed,
        enforce=cfg.enforce,
        calibration_enabled=cfg.calibration_enabled,
        init_design=cfg.init_design,
        full_confirm=cfg.full_confirm,
        solver_modes=cfg.solver_modes,
        dt_max_orbit_values=cfg.dt_max_orbit_values,
        dt_max_atm_values=cfg.dt_max_atm_values,
        maxiters_quick=cfg.maxiters_quick,
        maxiters_full=cfg.maxiters_full,
        scenario_weights=cfg.scenario_weights,
    )
end

@inline function _with_cfg(cfg::TunerConfig; kwargs...)
    return TunerConfig(; _tuner_config_fields(cfg)..., kwargs...)
end

Base.@kwdef struct ParameterSpec
    name::String
    kind::Symbol
    lower::Float64 = 0.0
    upper::Float64 = 1.0
    choices::Vector{String} = String[]
end

Base.@kwdef struct TuneCandidate
    id::Int
    values::Dict{String, Any}
    stage::String = "quick"
end

Base.@kwdef struct EvalTask
    candidate::TuneCandidate
    draw::Union{Nothing, Int} = nothing
end

Base.@kwdef struct SurrogateModel
    X::Matrix{Float64}
    y::Vector{Float64}
    Ainv::Matrix{Float64}
    alpha::Vector{Float64}
    length_scale::Float64
    residual_sigma::Float64
end

@inline function _parse_bool(raw::AbstractString)::Bool
    token = lowercase(strip(raw))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean token '$raw'."))
end

@inline function _parse_int(raw::AbstractString, label::String)::Int
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("Expected integer for $label, got '$raw'."))
    end
    return value
end

@inline function _parse_float(raw::AbstractString, label::String)::Float64
    value = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("Expected float for $label, got '$raw'."))
    end
    return value
end

@inline function _parse_float_list(raw::AbstractString, label::String)::Vector{Float64}
    vals = Float64[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(vals, _parse_float(t, label))
    end
    isempty(vals) && throw(ArgumentError("Expected non-empty float list for $label."))
    return vals
end

@inline function _parse_string_list(raw::AbstractString, label::String)::Vector{String}
    vals = String[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(vals, t)
    end
    isempty(vals) && throw(ArgumentError("Expected non-empty string list for $label."))
    return vals
end

function _parse_weights(raw::AbstractString)::Dict{String, Float64}
    out = Dict{String, Float64}()
    for token in split(raw, ",")
        pair = strip(token)
        isempty(pair) && continue
        parts = split(pair, ":", limit=2)
        length(parts) == 2 || throw(ArgumentError("Invalid scenario weight token '$pair'; expected scenario:weight."))
        name = strip(parts[1])
        val = _parse_float(strip(parts[2]), "scenario-weights")
        out[name] = val
    end
    isempty(out) && throw(ArgumentError("Scenario weights cannot be empty."))
    return out
end

function parse_cli(args::Vector{String})::TunerConfig
    cfg = TunerConfig()

    for arg in args
        if startswith(arg, "--base-manifest=")
            cfg = _with_cfg(cfg; base_manifest=abspath(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--outdir=")
            cfg = _with_cfg(cfg; outdir=abspath(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--n-init=")
            cfg = _with_cfg(cfg; n_init=_parse_int(split(arg, "=", limit=2)[2], "n-init"))
        elseif startswith(arg, "--n-bo-iters=")
            cfg = _with_cfg(cfg; n_bo_iters=_parse_int(split(arg, "=", limit=2)[2], "n-bo-iters"))
        elseif startswith(arg, "--batch-size=")
            cfg = _with_cfg(cfg; batch_size=_parse_int(split(arg, "=", limit=2)[2], "batch-size"))
        elseif startswith(arg, "--bo-pool-size=")
            cfg = _with_cfg(cfg; bo_pool_size=_parse_int(split(arg, "=", limit=2)[2], "bo-pool-size"))
        elseif startswith(arg, "--acq=")
            acq_raw = lowercase(strip(split(arg, "=", limit=2)[2]))
            acq_raw in ("lcb", "ei") || throw(ArgumentError("acq must be lcb|ei."))
            cfg = _with_cfg(cfg; acquisition=Symbol(acq_raw))
        elseif startswith(arg, "--kappa-start=")
            cfg = _with_cfg(cfg; kappa_start=_parse_float(split(arg, "=", limit=2)[2], "kappa-start"))
        elseif startswith(arg, "--kappa-end=")
            cfg = _with_cfg(cfg; kappa_end=_parse_float(split(arg, "=", limit=2)[2], "kappa-end"))
        elseif startswith(arg, "--ei-xi=")
            cfg = _with_cfg(cfg; ei_xi=_parse_float(split(arg, "=", limit=2)[2], "ei-xi"))
        elseif startswith(arg, "--surrogate-length-scale=")
            cfg = _with_cfg(cfg; surrogate_length_scale=_parse_float(split(arg, "=", limit=2)[2], "surrogate-length-scale"))
        elseif startswith(arg, "--surrogate-noise=")
            cfg = _with_cfg(cfg; surrogate_noise=_parse_float(split(arg, "=", limit=2)[2], "surrogate-noise"))
        elseif startswith(arg, "--local-topk=")
            cfg = _with_cfg(cfg; local_seed_topk=_parse_int(split(arg, "=", limit=2)[2], "local-topk"))
        elseif startswith(arg, "--local-steps=")
            cfg = _with_cfg(cfg; local_steps=_parse_int(split(arg, "=", limit=2)[2], "local-steps"))
        elseif startswith(arg, "--local-step-init=")
            cfg = _with_cfg(cfg; local_step_init=_parse_float(split(arg, "=", limit=2)[2], "local-step-init"))
        elseif startswith(arg, "--local-step-shrink=")
            cfg = _with_cfg(cfg; local_step_shrink=_parse_float(split(arg, "=", limit=2)[2], "local-step-shrink"))
        elseif startswith(arg, "--local-step-expand=")
            cfg = _with_cfg(cfg; local_step_expand=_parse_float(split(arg, "=", limit=2)[2], "local-step-expand"))
        elseif startswith(arg, "--local-min-improve=")
            cfg = _with_cfg(cfg; local_min_improve=_parse_float(split(arg, "=", limit=2)[2], "local-min-improve"))
        elseif startswith(arg, "--finalists=")
            cfg = _with_cfg(cfg; finalists=_parse_int(split(arg, "=", limit=2)[2], "finalists"))
        elseif startswith(arg, "--robust-draws=")
            cfg = _with_cfg(cfg; robust_draws=_parse_int(split(arg, "=", limit=2)[2], "robust-draws"))
        elseif startswith(arg, "--robust-alpha=")
            cfg = _with_cfg(cfg; robust_alpha=_parse_float(split(arg, "=", limit=2)[2], "robust-alpha"))
        elseif startswith(arg, "--robust-beta=")
            cfg = _with_cfg(cfg; robust_beta=_parse_float(split(arg, "=", limit=2)[2], "robust-beta"))
        elseif startswith(arg, "--uncertainty-atm-scale=")
            cfg = _with_cfg(cfg; uncertainty_atm_scale=_parse_float(split(arg, "=", limit=2)[2], "uncertainty-atm-scale"))
        elseif startswith(arg, "--uncertainty-ic-scale=")
            cfg = _with_cfg(cfg; uncertainty_ic_scale=_parse_float(split(arg, "=", limit=2)[2], "uncertainty-ic-scale"))
        elseif startswith(arg, "--telemetry-noise-frac=")
            cfg = _with_cfg(cfg; telemetry_noise_frac=_parse_float(split(arg, "=", limit=2)[2], "telemetry-noise-frac"))
        elseif startswith(arg, "--huber-delta=")
            cfg = _with_cfg(cfg; huber_delta=_parse_float(split(arg, "=", limit=2)[2], "huber-delta"))
        elseif startswith(arg, "--lambda-fail=")
            cfg = _with_cfg(cfg; lambda_fail=_parse_float(split(arg, "=", limit=2)[2], "lambda-fail"))
        elseif startswith(arg, "--lambda-time=")
            cfg = _with_cfg(cfg; lambda_time=_parse_float(split(arg, "=", limit=2)[2], "lambda-time"))
        elseif startswith(arg, "--runtime-budget-s=")
            cfg = _with_cfg(cfg; runtime_budget_s=_parse_float(split(arg, "=", limit=2)[2], "runtime-budget-s"))
        elseif startswith(arg, "--parallel-candidates=")
            cfg = _with_cfg(cfg; parallel_candidates=_parse_int(split(arg, "=", limit=2)[2], "parallel-candidates"))
        elseif startswith(arg, "--seed=")
            cfg = _with_cfg(cfg; seed=_parse_int(split(arg, "=", limit=2)[2], "seed"))
        elseif startswith(arg, "--enforce=")
            cfg = _with_cfg(cfg; enforce=_parse_bool(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--calibration=")
            cfg = _with_cfg(cfg; calibration_enabled=_parse_bool(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--init-design=")
            token = lowercase(strip(split(arg, "=", limit=2)[2]))
            token in ("lhs", "random") || throw(ArgumentError("init-design must be lhs|random."))
            cfg = _with_cfg(cfg; init_design=Symbol(token))
        elseif startswith(arg, "--full-confirm=")
            cfg = _with_cfg(cfg; full_confirm=_parse_bool(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--solver-modes=")
            cfg = _with_cfg(cfg; solver_modes=_parse_string_list(split(arg, "=", limit=2)[2], "solver-modes"))
        elseif startswith(arg, "--dt-orbit-values=")
            cfg = _with_cfg(cfg; dt_max_orbit_values=_parse_float_list(split(arg, "=", limit=2)[2], "dt-orbit-values"))
        elseif startswith(arg, "--dt-atm-values=")
            cfg = _with_cfg(cfg; dt_max_atm_values=_parse_float_list(split(arg, "=", limit=2)[2], "dt-atm-values"))
        elseif startswith(arg, "--maxiters-quick=")
            cfg = _with_cfg(cfg; maxiters_quick=_parse_int(split(arg, "=", limit=2)[2], "maxiters-quick"))
        elseif startswith(arg, "--maxiters-full=")
            cfg = _with_cfg(cfg; maxiters_full=_parse_int(split(arg, "=", limit=2)[2], "maxiters-full"))
        elseif startswith(arg, "--scenario-weights=")
            cfg = _with_cfg(cfg; scenario_weights=_parse_weights(split(arg, "=", limit=2)[2]))
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: --base-manifest, --outdir, --n-init, --n-bo-iters, --batch-size, " *
                "--bo-pool-size, --acq, --kappa-start, --kappa-end, --ei-xi, --surrogate-length-scale, --surrogate-noise, " *
                "--local-topk, --local-steps, --local-step-init, --local-step-shrink, --local-step-expand, --local-min-improve, " *
                "--finalists, --robust-draws, --robust-alpha, --robust-beta, --uncertainty-atm-scale, --uncertainty-ic-scale, " *
                "--telemetry-noise-frac, --huber-delta, --lambda-fail, --lambda-time, --runtime-budget-s, --parallel-candidates, " *
                "--seed, --enforce, --calibration, --init-design, --full-confirm, --solver-modes, --dt-orbit-values, --dt-atm-values, " *
                "--maxiters-quick, --maxiters-full, --scenario-weights"
            ))
        end
    end

    cfg.n_init > 0 || throw(ArgumentError("n-init must be > 0."))
    cfg.n_bo_iters >= 0 || throw(ArgumentError("n-bo-iters must be >= 0."))
    cfg.batch_size > 0 || throw(ArgumentError("batch-size must be > 0."))
    cfg.bo_pool_size > 0 || throw(ArgumentError("bo-pool-size must be > 0."))
    cfg.parallel_candidates > 0 || throw(ArgumentError("parallel-candidates must be > 0."))
    cfg.local_seed_topk > 0 || throw(ArgumentError("local-topk must be > 0."))
    cfg.local_steps >= 0 || throw(ArgumentError("local-steps must be >= 0."))
    cfg.local_step_init > 0.0 || throw(ArgumentError("local-step-init must be > 0."))
    cfg.local_step_shrink > 0.0 || throw(ArgumentError("local-step-shrink must be > 0."))
    cfg.local_step_expand > 0.0 || throw(ArgumentError("local-step-expand must be > 0."))
    cfg.finalists > 0 || throw(ArgumentError("finalists must be > 0."))
    cfg.robust_draws > 0 || throw(ArgumentError("robust-draws must be > 0."))
    cfg.runtime_budget_s > 0.0 || throw(ArgumentError("runtime-budget-s must be > 0."))
    all(v -> v > 0.0, cfg.dt_max_orbit_values) || throw(ArgumentError("dt-orbit-values must be > 0."))
    all(v -> v > 0.0, cfg.dt_max_atm_values) || throw(ArgumentError("dt-atm-values must be > 0."))
    return cfg
end

@inline function _progress_print(lock_ref::Union{Nothing, ReentrantLock}, msg::String)
    if lock_ref === nothing
        println(msg)
    else
        lock(lock_ref) do
            println(msg)
        end
    end
end

@inline function _env_with_default(name::String, default::String)::String
    raw = strip(get(ENV, name, ""))
    return isempty(raw) ? default : raw
end

function _candidate_runtime_policy_env_pairs()::Vector{Pair{String, String}}
    return Pair{String, String}[
        "OPENBLAS_NUM_THREADS" => _env_with_default("OPENBLAS_NUM_THREADS", "1"),
        "SPACEAGORA_INNER_THREAD_BUDGET" => _env_with_default("SPACEAGORA_INNER_THREAD_BUDGET", "1"),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => _env_with_default("SPACEAGORA_PARALLEL_POLICY_ADAPTIVE", "0"),
        "SPACEAGORA_EFFECTOR_PARALLEL" => _env_with_default("SPACEAGORA_EFFECTOR_PARALLEL", "off"),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => _env_with_default("SPACEAGORA_DENSITY_CALLBACK_PARALLEL", "off"),
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => _env_with_default("SPACEAGORA_CONTROL_CALLBACK_PARALLEL", "off"),
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => _env_with_default("SPACEAGORA_THERMAL_CALLBACK_PARALLEL", "off"),
        "SPACEAGORA_MULTIBODY_PARALLEL" => _env_with_default("SPACEAGORA_MULTIBODY_PARALLEL", "off"),
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => _env_with_default("SPACEAGORA_OUTER_PARALLEL_ACTIVE", "0")
    ]
end

@inline function _wrap_0_360(deg::Float64)::Float64
    out = mod(deg, 360.0)
    out < 0.0 && return out + 360.0
    return out
end

@inline _clamp_positive(x::Float64) = max(1e-9, x)
