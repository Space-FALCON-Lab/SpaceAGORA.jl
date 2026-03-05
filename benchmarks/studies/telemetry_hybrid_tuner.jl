const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_BASE_MANIFEST = joinpath(REPO_ROOT, "test", "telemetry_benchmark_manifest.toml")
const DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "telemetry_tuning", "hybrid")
const TELEMETRY_STUDY_SCRIPT = joinpath(REPO_ROOT, "src", "analysis", "verification", "telemetry_verification.jl")
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
            cfg = TunerConfig(; cfg..., base_manifest=abspath(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--outdir=")
            cfg = TunerConfig(; cfg..., outdir=abspath(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--n-init=")
            cfg = TunerConfig(; cfg..., n_init=_parse_int(split(arg, "=", limit=2)[2], "n-init"))
        elseif startswith(arg, "--n-bo-iters=")
            cfg = TunerConfig(; cfg..., n_bo_iters=_parse_int(split(arg, "=", limit=2)[2], "n-bo-iters"))
        elseif startswith(arg, "--batch-size=")
            cfg = TunerConfig(; cfg..., batch_size=_parse_int(split(arg, "=", limit=2)[2], "batch-size"))
        elseif startswith(arg, "--bo-pool-size=")
            cfg = TunerConfig(; cfg..., bo_pool_size=_parse_int(split(arg, "=", limit=2)[2], "bo-pool-size"))
        elseif startswith(arg, "--acq=")
            acq_raw = lowercase(strip(split(arg, "=", limit=2)[2]))
            acq_raw in ("lcb", "ei") || throw(ArgumentError("acq must be lcb|ei."))
            cfg = TunerConfig(; cfg..., acquisition=Symbol(acq_raw))
        elseif startswith(arg, "--kappa-start=")
            cfg = TunerConfig(; cfg..., kappa_start=_parse_float(split(arg, "=", limit=2)[2], "kappa-start"))
        elseif startswith(arg, "--kappa-end=")
            cfg = TunerConfig(; cfg..., kappa_end=_parse_float(split(arg, "=", limit=2)[2], "kappa-end"))
        elseif startswith(arg, "--ei-xi=")
            cfg = TunerConfig(; cfg..., ei_xi=_parse_float(split(arg, "=", limit=2)[2], "ei-xi"))
        elseif startswith(arg, "--surrogate-length-scale=")
            cfg = TunerConfig(; cfg..., surrogate_length_scale=_parse_float(split(arg, "=", limit=2)[2], "surrogate-length-scale"))
        elseif startswith(arg, "--surrogate-noise=")
            cfg = TunerConfig(; cfg..., surrogate_noise=_parse_float(split(arg, "=", limit=2)[2], "surrogate-noise"))
        elseif startswith(arg, "--local-topk=")
            cfg = TunerConfig(; cfg..., local_seed_topk=_parse_int(split(arg, "=", limit=2)[2], "local-topk"))
        elseif startswith(arg, "--local-steps=")
            cfg = TunerConfig(; cfg..., local_steps=_parse_int(split(arg, "=", limit=2)[2], "local-steps"))
        elseif startswith(arg, "--local-step-init=")
            cfg = TunerConfig(; cfg..., local_step_init=_parse_float(split(arg, "=", limit=2)[2], "local-step-init"))
        elseif startswith(arg, "--local-step-shrink=")
            cfg = TunerConfig(; cfg..., local_step_shrink=_parse_float(split(arg, "=", limit=2)[2], "local-step-shrink"))
        elseif startswith(arg, "--local-step-expand=")
            cfg = TunerConfig(; cfg..., local_step_expand=_parse_float(split(arg, "=", limit=2)[2], "local-step-expand"))
        elseif startswith(arg, "--local-min-improve=")
            cfg = TunerConfig(; cfg..., local_min_improve=_parse_float(split(arg, "=", limit=2)[2], "local-min-improve"))
        elseif startswith(arg, "--finalists=")
            cfg = TunerConfig(; cfg..., finalists=_parse_int(split(arg, "=", limit=2)[2], "finalists"))
        elseif startswith(arg, "--robust-draws=")
            cfg = TunerConfig(; cfg..., robust_draws=_parse_int(split(arg, "=", limit=2)[2], "robust-draws"))
        elseif startswith(arg, "--robust-alpha=")
            cfg = TunerConfig(; cfg..., robust_alpha=_parse_float(split(arg, "=", limit=2)[2], "robust-alpha"))
        elseif startswith(arg, "--robust-beta=")
            cfg = TunerConfig(; cfg..., robust_beta=_parse_float(split(arg, "=", limit=2)[2], "robust-beta"))
        elseif startswith(arg, "--uncertainty-atm-scale=")
            cfg = TunerConfig(; cfg..., uncertainty_atm_scale=_parse_float(split(arg, "=", limit=2)[2], "uncertainty-atm-scale"))
        elseif startswith(arg, "--uncertainty-ic-scale=")
            cfg = TunerConfig(; cfg..., uncertainty_ic_scale=_parse_float(split(arg, "=", limit=2)[2], "uncertainty-ic-scale"))
        elseif startswith(arg, "--telemetry-noise-frac=")
            cfg = TunerConfig(; cfg..., telemetry_noise_frac=_parse_float(split(arg, "=", limit=2)[2], "telemetry-noise-frac"))
        elseif startswith(arg, "--huber-delta=")
            cfg = TunerConfig(; cfg..., huber_delta=_parse_float(split(arg, "=", limit=2)[2], "huber-delta"))
        elseif startswith(arg, "--lambda-fail=")
            cfg = TunerConfig(; cfg..., lambda_fail=_parse_float(split(arg, "=", limit=2)[2], "lambda-fail"))
        elseif startswith(arg, "--lambda-time=")
            cfg = TunerConfig(; cfg..., lambda_time=_parse_float(split(arg, "=", limit=2)[2], "lambda-time"))
        elseif startswith(arg, "--runtime-budget-s=")
            cfg = TunerConfig(; cfg..., runtime_budget_s=_parse_float(split(arg, "=", limit=2)[2], "runtime-budget-s"))
        elseif startswith(arg, "--parallel-candidates=")
            cfg = TunerConfig(; cfg..., parallel_candidates=_parse_int(split(arg, "=", limit=2)[2], "parallel-candidates"))
        elseif startswith(arg, "--seed=")
            cfg = TunerConfig(; cfg..., seed=_parse_int(split(arg, "=", limit=2)[2], "seed"))
        elseif startswith(arg, "--enforce=")
            cfg = TunerConfig(; cfg..., enforce=_parse_bool(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--calibration=")
            cfg = TunerConfig(; cfg..., calibration_enabled=_parse_bool(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--init-design=")
            token = lowercase(strip(split(arg, "=", limit=2)[2]))
            token in ("lhs", "random") || throw(ArgumentError("init-design must be lhs|random."))
            cfg = TunerConfig(; cfg..., init_design=Symbol(token))
        elseif startswith(arg, "--full-confirm=")
            cfg = TunerConfig(; cfg..., full_confirm=_parse_bool(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--solver-modes=")
            cfg = TunerConfig(; cfg..., solver_modes=_parse_string_list(split(arg, "=", limit=2)[2], "solver-modes"))
        elseif startswith(arg, "--dt-orbit-values=")
            cfg = TunerConfig(; cfg..., dt_max_orbit_values=_parse_float_list(split(arg, "=", limit=2)[2], "dt-orbit-values"))
        elseif startswith(arg, "--dt-atm-values=")
            cfg = TunerConfig(; cfg..., dt_max_atm_values=_parse_float_list(split(arg, "=", limit=2)[2], "dt-atm-values"))
        elseif startswith(arg, "--maxiters-quick=")
            cfg = TunerConfig(; cfg..., maxiters_quick=_parse_int(split(arg, "=", limit=2)[2], "maxiters-quick"))
        elseif startswith(arg, "--maxiters-full=")
            cfg = TunerConfig(; cfg..., maxiters_full=_parse_int(split(arg, "=", limit=2)[2], "maxiters-full"))
        elseif startswith(arg, "--scenario-weights=")
            cfg = TunerConfig(; cfg..., scenario_weights=_parse_weights(split(arg, "=", limit=2)[2]))
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

function _parameter_space(cfg::TunerConfig)::Vector{ParameterSpec}
    return ParameterSpec[
        ParameterSpec(name="epoch_shift_s", kind=:integer, lower=-180.0, upper=180.0),
        ParameterSpec(name="ra_scale", kind=:continuous, lower=0.995, upper=1.005),
        ParameterSpec(name="rp_altitude_offset_m", kind=:continuous, lower=-3000.0, upper=3000.0),
        ParameterSpec(name="i_offset_deg", kind=:continuous, lower=-0.03, upper=0.03),
        ParameterSpec(name="aop_offset_deg", kind=:continuous, lower=-0.2, upper=0.2),
        ParameterSpec(name="raan_offset_deg", kind=:continuous, lower=-0.2, upper=0.2),
        ParameterSpec(name="ta_offset_deg", kind=:continuous, lower=-0.5, upper=0.5),
        ParameterSpec(name="bus_mass_scale", kind=:continuous, lower=0.99, upper=1.01),
        ParameterSpec(name="prop_mass_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="panel_mass_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="bus_dims_scale", kind=:continuous, lower=0.995, upper=1.005),
        ParameterSpec(name="panel_dims_scale", kind=:continuous, lower=0.995, upper=1.005),
        ParameterSpec(name="panel_offset_scale", kind=:continuous, lower=0.99, upper=1.01),
        ParameterSpec(name="srp_cr_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="srp_area_scale", kind=:continuous, lower=0.95, upper=1.05),
        ParameterSpec(name="gravity_variant", kind=:categorical, choices=copy(GRAVITY_VARIANTS)),
        ParameterSpec(name="dv_global_scale", kind=:continuous, lower=0.8, upper=1.2),
        ParameterSpec(name="dv_early_scale", kind=:continuous, lower=0.8, upper=1.2),
        ParameterSpec(name="dv_orbit7_bias_mps", kind=:continuous, lower=-0.3, upper=0.3),
        ParameterSpec(name="solver_mode", kind=:categorical, choices=copy(cfg.solver_modes)),
        ParameterSpec(name="dt_max_orbit_s", kind=:categorical, choices=string.(cfg.dt_max_orbit_values)),
        ParameterSpec(name="dt_max_atm_s", kind=:categorical, choices=string.(cfg.dt_max_atm_values))
    ]
end

@inline function _sample_value(rng::AbstractRNG, p::ParameterSpec)
    if p.kind == :continuous
        return rand(rng) * (p.upper - p.lower) + p.lower
    elseif p.kind == :integer
        lo = ceil(Int, p.lower)
        hi = floor(Int, p.upper)
        lo <= hi || throw(ArgumentError("integer parameter $(p.name) has empty range."))
        return rand(rng, lo:hi)
    elseif p.kind == :categorical
        return rand(rng, p.choices)
    end
    throw(ArgumentError("Unsupported parameter kind $(p.kind)."))
end

function _candidate_signature(c::TuneCandidate)::String
    chunks = String[]
    for key in PARAMETER_ORDER
        push!(chunks, string(key, "=", c.values[key]))
    end
    return join(chunks, "|")
end

function _candidate_row_nt(c::TuneCandidate)
    return (
        candidate_id=c.id,
        epoch_shift_s=Int(c.values["epoch_shift_s"]),
        ra_scale=Float64(c.values["ra_scale"]),
        rp_altitude_offset_m=Float64(c.values["rp_altitude_offset_m"]),
        i_offset_deg=Float64(c.values["i_offset_deg"]),
        aop_offset_deg=Float64(c.values["aop_offset_deg"]),
        raan_offset_deg=Float64(c.values["raan_offset_deg"]),
        ta_offset_deg=Float64(c.values["ta_offset_deg"]),
        bus_mass_scale=Float64(c.values["bus_mass_scale"]),
        prop_mass_scale=Float64(c.values["prop_mass_scale"]),
        panel_mass_scale=Float64(c.values["panel_mass_scale"]),
        bus_dims_scale=Float64(c.values["bus_dims_scale"]),
        panel_dims_scale=Float64(c.values["panel_dims_scale"]),
        panel_offset_scale=Float64(c.values["panel_offset_scale"]),
        srp_cr_scale=Float64(c.values["srp_cr_scale"]),
        srp_area_scale=Float64(c.values["srp_area_scale"]),
        gravity_variant=String(c.values["gravity_variant"]),
        dv_global_scale=Float64(c.values["dv_global_scale"]),
        dv_early_scale=Float64(c.values["dv_early_scale"]),
        dv_orbit7_bias_mps=Float64(c.values["dv_orbit7_bias_mps"]),
        solver_mode_requested=String(c.values["solver_mode"]),
        dt_max_orbit_requested_s=Float64(c.values["dt_max_orbit_s"]),
        dt_max_atm_requested_s=Float64(c.values["dt_max_atm_s"])
    )
end

@inline function _candidate_vector(params::Vector{ParameterSpec}, c::TuneCandidate)::Vector{Float64}
    vec = Vector{Float64}(undef, length(params))
    for i in eachindex(params)
        p = params[i]
        v = c.values[p.name]
        if p.kind == :continuous || p.kind == :integer
            span = p.upper - p.lower
            vec[i] = span > 0 ? clamp((Float64(v) - p.lower) / span, 0.0, 1.0) : 0.5
        else
            idx = findfirst(==(String(v)), p.choices)
            if idx === nothing || length(p.choices) <= 1
                vec[i] = 0.0
            else
                vec[i] = (idx - 1) / (length(p.choices) - 1)
            end
        end
    end
    return vec
end

function _candidate_from_vector(
    params::Vector{ParameterSpec},
    u::Vector{Float64},
    id::Int,
    stage::String
)::TuneCandidate
    values = Dict{String, Any}()
    for i in eachindex(params)
        p = params[i]
        ui = clamp(u[i], 0.0, 1.0)
        if p.kind == :continuous
            values[p.name] = p.lower + ui * (p.upper - p.lower)
        elseif p.kind == :integer
            lo = ceil(Int, p.lower)
            hi = floor(Int, p.upper)
            x = lo + floor(Int, ui * (hi - lo + 1))
            values[p.name] = clamp(x, lo, hi)
        else
            idx = clamp(1 + floor(Int, ui * length(p.choices)), 1, length(p.choices))
            values[p.name] = p.choices[idx]
        end
    end

    values["dt_max_orbit_s"] = Float64(values["dt_max_orbit_s"])
    values["dt_max_atm_s"] = Float64(values["dt_max_atm_s"])
    values["epoch_shift_s"] = Int(values["epoch_shift_s"])
    return TuneCandidate(id=id, values=values, stage=stage)
end

function _lhs_unit_matrix(rng::AbstractRNG, n::Int, d::Int)::Matrix{Float64}
    m = Matrix{Float64}(undef, n, d)
    for j in 1:d
        perm = randperm(rng, n)
        for i in 1:n
            m[i, j] = (perm[i] - rand(rng)) / n
        end
    end
    return m
end

function _initial_design(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    n::Int;
    start_id::Int=1,
    stage::String="quick_global_init"
)::Vector{TuneCandidate}
    rng = MersenneTwister(hash((cfg.seed, stage, "initial_design")))
    out = TuneCandidate[]
    seen = Set{String}()

    if cfg.init_design == :lhs
        u = _lhs_unit_matrix(rng, n, length(params))
        next_id = start_id
        for i in 1:n
            cand = _candidate_from_vector(params, vec(u[i, :]), next_id, stage)
            sig = _candidate_signature(cand)
            if sig in seen
                continue
            end
            push!(seen, sig)
            push!(out, cand)
            next_id += 1
        end
        while length(out) < n
            vals = Dict{String, Any}()
            for p in params
                vals[p.name] = _sample_value(rng, p)
            end
            vals["dt_max_orbit_s"] = Float64(vals["dt_max_orbit_s"])
            vals["dt_max_atm_s"] = Float64(vals["dt_max_atm_s"])
            vals["epoch_shift_s"] = Int(vals["epoch_shift_s"])
            cand = TuneCandidate(id=start_id + length(out), values=vals, stage=stage)
            sig = _candidate_signature(cand)
            if sig in seen
                continue
            end
            push!(seen, sig)
            push!(out, cand)
        end
        return out
    end

    next_id = start_id
    while length(out) < n
        vals = Dict{String, Any}()
        for p in params
            vals[p.name] = _sample_value(rng, p)
        end
        vals["dt_max_orbit_s"] = Float64(vals["dt_max_orbit_s"])
        vals["dt_max_atm_s"] = Float64(vals["dt_max_atm_s"])
        vals["epoch_shift_s"] = Int(vals["epoch_shift_s"])
        cand = TuneCandidate(id=next_id, values=vals, stage=stage)
        sig = _candidate_signature(cand)
        if sig in seen
            next_id += 1
            continue
        end
        push!(seen, sig)
        push!(out, cand)
        next_id += 1
    end
    return out
end

@inline function _rbf_kernel(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, l::Float64)::Float64
    d2 = 0.0
    @inbounds for i in eachindex(x)
        d = x[i] - y[i]
        d2 += d * d
    end
    return exp(-d2 / (2.0 * l * l))
end

@inline function _phi(z::Float64)::Float64
    return exp(-0.5 * z * z) / sqrt(2.0 * π)
end

@inline function _Phi(z::Float64)::Float64
    return 0.5 * (1.0 + tanh(sqrt(2.0 / π) * (z + 0.044715 * z^3)))
end

function _fit_surrogate(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    candidates::Vector{TuneCandidate},
    y::Vector{Float64}
)::Union{Nothing, SurrogateModel}
    n = length(candidates)
    n == 0 && return nothing

    X = Matrix{Float64}(undef, n, length(params))
    for i in 1:n
        X[i, :] = _candidate_vector(params, candidates[i])
    end

    l = cfg.surrogate_length_scale
    noise = max(cfg.surrogate_noise, 1.0e-9)

    K = Matrix{Float64}(undef, n, n)
    for i in 1:n
        xi = @view X[i, :]
        for j in i:n
            xj = @view X[j, :]
            kval = _rbf_kernel(xi, xj, l)
            K[i, j] = kval
            K[j, i] = kval
        end
    end

    A = K + noise * I
    Ainv = Matrix(inv(Matrix(A)))
    alpha = Ainv * y

    diagA = diag(Ainv)
    loo_pred = similar(y)
    for i in eachindex(y)
        if abs(diagA[i]) < 1.0e-12
            loo_pred[i] = y[i]
        else
            loo_pred[i] = y[i] - alpha[i] / diagA[i]
        end
    end
    residual_sigma = length(y) > 2 ? std(y .- loo_pred) : (length(y) == 1 ? 0.0 : std(y))

    return SurrogateModel(
        X=X,
        y=y,
        Ainv=Ainv,
        alpha=alpha,
        length_scale=l,
        residual_sigma=max(residual_sigma, 1.0e-9)
    )
end

function _predict(model::SurrogateModel, x::Vector{Float64})::Tuple{Float64, Float64}
    n = size(model.X, 1)
    k = Vector{Float64}(undef, n)
    dmin2 = Inf
    for i in 1:n
        xi = @view model.X[i, :]
        k[i] = _rbf_kernel(x, xi, model.length_scale)
        d2 = sum((x .- xi) .^ 2)
        d2 < dmin2 && (dmin2 = d2)
    end

    mu = dot(k, model.alpha)
    var = 1.0 - dot(k, model.Ainv * k)
    base_sigma = sqrt(max(var, 1.0e-12))
    density = exp(-dmin2 / (2.0 * model.length_scale * model.length_scale))
    sigma = base_sigma + model.residual_sigma * (1.0 - density)
    return mu, max(sigma, 1.0e-12)
end

@inline function _ei(mu::Float64, sigma::Float64, best::Float64, xi::Float64)::Float64
    sigma <= 1.0e-12 && return 0.0
    improve = best - mu - xi
    z = improve / sigma
    return improve * _Phi(z) + sigma * _phi(z)
end

@inline function _kappa_schedule(cfg::TunerConfig, iter::Int)::Float64
    cfg.n_bo_iters <= 1 && return cfg.kappa_end
    t = (iter - 1) / max(cfg.n_bo_iters - 1, 1)
    return cfg.kappa_start + (cfg.kappa_end - cfg.kappa_start) * t
end

function _propose_batch(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    model::SurrogateModel,
    observed_scores::Vector{Float64},
    seen::Set{String},
    next_id::Int,
    iter::Int,
    q::Int
)::Vector{TuneCandidate}
    proposed = TuneCandidate[]
    local_seen = Set{String}()
    best_obs = minimum(observed_scores)
    kappa = _kappa_schedule(cfg, iter)

    for j in 1:q
        rng = MersenneTwister(hash((cfg.seed, "bo_pool", iter, j, next_id)))
        pool = TuneCandidate[]
        attempts = 0
        max_attempts = max(cfg.bo_pool_size * 40, 1024)
        while length(pool) < cfg.bo_pool_size && attempts < max_attempts
            vals = Dict{String, Any}()
            for p in params
                vals[p.name] = _sample_value(rng, p)
            end
            vals["dt_max_orbit_s"] = Float64(vals["dt_max_orbit_s"])
            vals["dt_max_atm_s"] = Float64(vals["dt_max_atm_s"])
            vals["epoch_shift_s"] = Int(vals["epoch_shift_s"])
            cand = TuneCandidate(id=next_id + j - 1, values=vals, stage="quick_global_bo")
            sig = _candidate_signature(cand)
            if sig in seen || sig in local_seen
                attempts += 1
                continue
            end
            push!(pool, cand)
            push!(local_seen, sig)
            attempts += 1
        end
        isempty(pool) && break

        best_idx = 1
        best_val = cfg.acquisition == :ei ? -Inf : Inf
        best_sig = _candidate_signature(pool[1])

        for i in eachindex(pool)
            cand = pool[i]
            x = _candidate_vector(params, cand)
            mu, sigma = _predict(model, x)
            acq_val = if cfg.acquisition == :ei
                _ei(mu, sigma, best_obs, cfg.ei_xi)
            else
                mu - kappa * sigma
            end

            sig = _candidate_signature(cand)
            if cfg.acquisition == :ei
                if acq_val > best_val || (acq_val == best_val && sig < best_sig)
                    best_idx = i
                    best_val = acq_val
                    best_sig = sig
                end
            else
                if acq_val < best_val || (acq_val == best_val && sig < best_sig)
                    best_idx = i
                    best_val = acq_val
                    best_sig = sig
                end
            end
        end

        chosen = pool[best_idx]
        push!(proposed, chosen)
        push!(seen, _candidate_signature(chosen))
    end

    return proposed
end

function _candidate_paths(outdir::String, stage::String, id::Int; draw::Union{Nothing, Int}=nothing)
    tag = draw === nothing ? @sprintf("cand_%04d_%s", id, stage) : @sprintf("cand_%04d_%s_draw%03d", id, stage, draw)
    manifests_dir = joinpath(outdir, "manifests")
    results_dir = joinpath(outdir, "results")
    logs_dir = joinpath(outdir, "logs")
    mkpath(manifests_dir)
    mkpath(results_dir)
    mkpath(logs_dir)
    return (
        manifest=joinpath(manifests_dir, "$(tag).toml"),
        summary=joinpath(results_dir, "$(tag)_summary.csv"),
        errors=joinpath(results_dir, "$(tag)_errors.csv"),
        log=joinpath(logs_dir, "$(tag).log")
    )
end

function _apply_uncertainty_to_scenario!(
    sc::Dict{String, Any},
    rng::AbstractRNG,
    cfg::TunerConfig
)
    if haskey(sc, "atmosphere_truth")
        at = Dict{String, Any}(sc["atmosphere_truth"])
        base_seed = Int(get(at, "gram_seed", 1001))
        at["gram_seed"] = base_seed + rand(rng, 1:1_000_000)

        base_scales = if haskey(at, "gram_perturbation_scales")
            [Float64(v) for v in at["gram_perturbation_scales"]]
        else
            [0.0, 0.0, 0.0, 0.0]
        end
        jitter = cfg.uncertainty_atm_scale
        pert = Float64[]
        for s in base_scales
            sigma = s > 0.0 ? jitter * s : jitter
            push!(pert, max(0.0, s + sigma * randn(rng)))
        end
        at["gram_perturbation_scales"] = pert
        sc["atmosphere_truth"] = at
    end

    frac = cfg.uncertainty_ic_scale
    if haskey(sc, "ra_m")
        sc["ra_m"] = _clamp_positive(Float64(sc["ra_m"]) * (1.0 + frac * randn(rng)))
    end
    if haskey(sc, "rp_altitude_m")
        sc["rp_altitude_m"] = _clamp_positive(Float64(sc["rp_altitude_m"]) + 1500.0 * frac * randn(rng))
    end
    if haskey(sc, "i_deg")
        sc["i_deg"] = clamp(Float64(sc["i_deg"]) + 0.5 * frac * randn(rng), 0.0, 180.0)
    end
    if haskey(sc, "aop_deg")
        sc["aop_deg"] = _wrap_0_360(Float64(sc["aop_deg"]) + 2.0 * frac * randn(rng))
    end
    if haskey(sc, "raan_deg")
        sc["raan_deg"] = _wrap_0_360(Float64(sc["raan_deg"]) + 2.0 * frac * randn(rng))
    end
    if haskey(sc, "ta_deg")
        sc["ta_deg"] = _wrap_0_360(Float64(sc["ta_deg"]) + 5.0 * frac * randn(rng))
    end

    return sc
end

function _apply_candidate_to_manifest(
    base_manifest::Dict{String, Any},
    cand::TuneCandidate,
    cfg::TunerConfig;
    uncertainty_draw::Union{Nothing, Int}=nothing
)::Dict{String, Any}
    doc = deepcopy(base_manifest)
    haskey(doc, "scenarios") || throw(ArgumentError("Manifest missing top-level scenarios array."))
    scenarios = doc["scenarios"]
    scenarios isa AbstractVector || throw(ArgumentError("Manifest scenarios must be an array."))

    shift_s = Int(cand.values["epoch_shift_s"])

    for i in eachindex(scenarios)
        sc = Dict{String, Any}(scenarios[i])

        if haskey(sc, "initial_time")
            init = Dict{String, Any}(sc["initial_time"])
            sec_raw = Float64(get(init, "second", 0.0))
            sec_int = floor(Int, sec_raw)
            sec_frac = sec_raw - sec_int
            sec_ms = round(Int, sec_frac * 1000)
            dt0 = DateTime(
                Int(init["year"]),
                Int(init["month"]),
                Int(init["day"]),
                Int(init["hour"]),
                Int(init["minute"]),
                sec_int,
                sec_ms
            )
            dt_shifted = dt0 + Second(shift_s)
            init["year"] = year(dt_shifted)
            init["month"] = month(dt_shifted)
            init["day"] = day(dt_shifted)
            init["hour"] = hour(dt_shifted)
            init["minute"] = minute(dt_shifted)
            init["second"] = second(dt_shifted) + millisecond(dt_shifted) / 1000
            sc["initial_time"] = init
        end

        if haskey(sc, "ra_m")
            sc["ra_m"] = _clamp_positive(Float64(sc["ra_m"]) * Float64(cand.values["ra_scale"]))
        end
        if haskey(sc, "rp_altitude_m")
            sc["rp_altitude_m"] = _clamp_positive(Float64(sc["rp_altitude_m"]) + Float64(cand.values["rp_altitude_offset_m"]))
        end
        if haskey(sc, "i_deg")
            sc["i_deg"] = clamp(Float64(sc["i_deg"]) + Float64(cand.values["i_offset_deg"]), 0.0, 180.0)
        end
        if haskey(sc, "aop_deg")
            sc["aop_deg"] = _wrap_0_360(Float64(sc["aop_deg"]) + Float64(cand.values["aop_offset_deg"]))
        end
        if haskey(sc, "raan_deg")
            sc["raan_deg"] = _wrap_0_360(Float64(sc["raan_deg"]) + Float64(cand.values["raan_offset_deg"]))
        end
        if haskey(sc, "ta_deg")
            sc["ta_deg"] = _wrap_0_360(Float64(sc["ta_deg"]) + Float64(cand.values["ta_offset_deg"]))
        end

        if haskey(sc, "spacecraft")
            spacecraft = Dict{String, Any}(sc["spacecraft"])
            haskey(spacecraft, "bus_mass_kg") && (spacecraft["bus_mass_kg"] = _clamp_positive(Float64(spacecraft["bus_mass_kg"]) * Float64(cand.values["bus_mass_scale"])))
            haskey(spacecraft, "prop_mass_kg") && (spacecraft["prop_mass_kg"] = _clamp_positive(Float64(spacecraft["prop_mass_kg"]) * Float64(cand.values["prop_mass_scale"])))
            haskey(spacecraft, "panel_mass_each_kg") && (spacecraft["panel_mass_each_kg"] = _clamp_positive(Float64(spacecraft["panel_mass_each_kg"]) * Float64(cand.values["panel_mass_scale"])))
            haskey(spacecraft, "panel_offset_y_m") && (spacecraft["panel_offset_y_m"] = _clamp_positive(Float64(spacecraft["panel_offset_y_m"]) * Float64(cand.values["panel_offset_scale"])))

            if haskey(spacecraft, "bus_dims_m")
                spacecraft["bus_dims_m"] = [Float64(v) * Float64(cand.values["bus_dims_scale"]) for v in spacecraft["bus_dims_m"]]
            end
            if haskey(spacecraft, "panel_dims_m")
                spacecraft["panel_dims_m"] = [Float64(v) * Float64(cand.values["panel_dims_scale"]) for v in spacecraft["panel_dims_m"]]
            end
            sc["spacecraft"] = spacecraft
        end

        haskey(sc, "srp_cr") && (sc["srp_cr"] = _clamp_positive(Float64(sc["srp_cr"]) * Float64(cand.values["srp_cr_scale"])))
        haskey(sc, "srp_area_m2") && (sc["srp_area_m2"] = _clamp_positive(Float64(sc["srp_area_m2"]) * Float64(cand.values["srp_area_scale"])))

        gravity_variant = String(cand.values["gravity_variant"])
        if gravity_variant == "j2_only"
            sc["gravity_model"] = "inverse_squared_j2"
            sc["gravity_harmonics_degree"] = 0
            sc["gravity_harmonics_order"] = 0
        else
            degree = if gravity_variant == "harm_deg2"
                2
            elseif gravity_variant == "harm_deg4"
                4
            elseif gravity_variant == "harm_deg8"
                8
            else
                20
            end
            sc["gravity_model"] = "inverse_squared"
            sc["gravity_harmonics_degree"] = degree
            sc["gravity_harmonics_order"] = degree
        end

        if haskey(sc, "maneuvers")
            mv = Dict{String, Any}(sc["maneuvers"])
            if haskey(mv, "orbit_numbers") && haskey(mv, "delta_v_mps")
                orbits = [Int(v) for v in mv["orbit_numbers"]]
                dvs = [Float64(v) for v in mv["delta_v_mps"]]
                tuned = similar(dvs)
                for j in eachindex(dvs)
                    scale = Float64(cand.values["dv_global_scale"])
                    if orbits[j] <= 50
                        scale *= Float64(cand.values["dv_early_scale"])
                    end
                    dv = dvs[j] * scale
                    if orbits[j] == 7
                        dv += Float64(cand.values["dv_orbit7_bias_mps"])
                    end
                    tuned[j] = dv
                end
                mv["delta_v_mps"] = tuned
                sc["maneuvers"] = mv
            end
        end

        if haskey(sc, "calibration")
            cal = Dict{String, Any}(sc["calibration"])
            cal["enabled"] = cfg.calibration_enabled
            sc["calibration"] = cal
        end

        if uncertainty_draw !== nothing
            rng = MersenneTwister(hash((cfg.seed, "uncertainty", String(get(sc, "name", "scenario")), cand.id, uncertainty_draw)))
            _apply_uncertainty_to_scenario!(sc, rng, cfg)
        end

        scenarios[i] = sc
    end

    doc["scenarios"] = scenarios
    return doc
end

@inline function _event_metric(df::DataFrame, scenario::String, event::String, col::Symbol)
    mask = (df.scenario .== scenario) .& (df.event .== event)
    idx = findfirst(mask)
    if idx === nothing
        return missing
    end
    value = df[idx, col]
    if value isa Missing
        return missing
    end
    return Float64(value)
end

@inline function _huber(x::Float64, delta::Float64)::Float64
    ax = abs(x)
    if ax <= delta
        return 0.5 * x * x
    end
    return delta * (ax - 0.5 * delta)
end

function _objective_from_summary(
    summary_df::DataFrame,
    cfg::TunerConfig;
    run_failed::Bool,
    runtime_s::Float64,
    noise_rng::Union{Nothing, AbstractRNG}=nothing
)
    if run_failed || nrow(summary_df) == 0
        rt_pen = cfg.lambda_time * max(0.0, runtime_s / cfg.runtime_budget_s - 1.0)
        score = cfg.lambda_fail + rt_pen
        return (
            objective=score,
            base_loss=0.0,
            fail_penalty=cfg.lambda_fail,
            runtime_penalty=rt_pen,
            all_pass=false,
            failed_rows=0,
            per_scenario=""
        )
    end

    base_loss = 0.0
    all_pass = true
    failed_rows = 0
    scenario_acc = Dict{String, Float64}()

    for row in eachrow(summary_df)
        scenario = String(row.scenario)
        w = get(cfg.scenario_weights, scenario, 1.0)

        rmse = Float64(row.rmse_km)
        absv = Float64(row.max_abs_km)
        lim_rmse = Float64(row.limit_max_rmse_km)
        lim_abs = Float64(row.limit_max_abs_km)

        rmse_ratio = if isfinite(lim_rmse) && lim_rmse > 0.0
            rmse / lim_rmse
        elseif isfinite(Float64(row.limit_nmae)) && Float64(row.limit_nmae) > 0.0
            Float64(row.nmae) / Float64(row.limit_nmae)
        else
            rmse
        end
        abs_ratio = absv / max(lim_abs, 1.0e-12)

        if noise_rng !== nothing
            rmse_ratio = max(0.0, rmse_ratio * (1.0 + cfg.telemetry_noise_frac * randn(noise_rng)))
            abs_ratio = max(0.0, abs_ratio * (1.0 + cfg.telemetry_noise_frac * randn(noise_rng)))
        end

        term = w * (_huber(rmse_ratio, cfg.huber_delta) + 0.5 * _huber(abs_ratio, cfg.huber_delta))
        base_loss += term
        scenario_acc[scenario] = get(scenario_acc, scenario, 0.0) + term

        pass_row = Bool(row.pass)
        all_pass &= pass_row
        pass_row || (failed_rows += 1)
    end

    fail_pen = cfg.lambda_fail * (run_failed ? 1.0 : 0.0)
    runtime_pen = cfg.lambda_time * max(0.0, runtime_s / cfg.runtime_budget_s - 1.0)
    score = base_loss + fail_pen + runtime_pen

    parts = [string(k, "=", @sprintf("%.6f", scenario_acc[k])) for k in sort(collect(keys(scenario_acc)))]
    return (
        objective=score,
        base_loss=base_loss,
        fail_penalty=fail_pen,
        runtime_penalty=runtime_pen,
        all_pass=all_pass,
        failed_rows=failed_rows,
        per_scenario=join(parts, ";")
    )
end

function _parse_run_metrics(summary_df::DataFrame, wall_s::Float64)
    if nrow(summary_df) == 0
        return (runtime_s=wall_s, scenario_count=0, event_count=0)
    end

    sim_runtime = try
        Float64(summary_df.total_runtime_s[1])
    catch
        wall_s
    end
    scenarios = length(unique(String.(summary_df.scenario)))
    events = nrow(summary_df)
    return (runtime_s=sim_runtime, scenario_count=scenarios, event_count=events)
end

function evaluate_candidate(
    cfg::TunerConfig,
    stage::String,
    profile::Symbol,
    cand::TuneCandidate,
    base_manifest::Dict{String, Any};
    uncertainty_draw::Union{Nothing, Int}=nothing,
    progress_lock::Union{Nothing, ReentrantLock}=nothing
)
    paths = _candidate_paths(cfg.outdir, stage, cand.id; draw=uncertainty_draw)
    manifest_doc = _apply_candidate_to_manifest(base_manifest, cand, cfg; uncertainty_draw=uncertainty_draw)
    open(paths.manifest, "w") do io
        TOML.print(io, manifest_doc)
    end

    cmd = `$(Base.julia_cmd()) --project=$(joinpath(REPO_ROOT, ".AGORA")) --startup-file=no $(TELEMETRY_STUDY_SCRIPT) --profile=$(String(profile)) --manifest=$(paths.manifest) --enforce=$(cfg.enforce ? "1" : "0") --out-summary=$(paths.summary) --out-errors=$(paths.errors)`

    env_pairs = Pair{String, String}[
        "SPACEAGORA_TELEMETRY_SOLVER_MODE" => String(cand.values["solver_mode"]),
        "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => string(Float64(cand.values["dt_max_orbit_s"])),
        "SPACEAGORA_TELEMETRY_DT_MAX_ATM" => string(Float64(cand.values["dt_max_atm_s"]))
    ]
    append!(env_pairs, _candidate_runtime_policy_env_pairs())
    if profile == :quick && cfg.maxiters_quick !== nothing
        push!(env_pairs, "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => string(cfg.maxiters_quick))
    elseif profile == :full && cfg.maxiters_full !== nothing
        push!(env_pairs, "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => string(cfg.maxiters_full))
    end

    run_ok = true
    error_text = ""
    wall_s = @elapsed begin
        open(paths.log, "w") do io
            try
                withenv(env_pairs...) do
                    proc = run(pipeline(cmd, stdout=io, stderr=io); wait=false)
                    started = time()
                    next_heartbeat = started + TUNER_HEARTBEAT_SECONDS
                    while process_running(proc)
                        sleep(2.0)
                        now_s = time()
                        if now_s >= next_heartbeat
                            elapsed_min = (now_s - started) / 60.0
                            _progress_print(
                                progress_lock,
                                @sprintf(
                                    "[heartbeat %s/%s] candidate=%d elapsed=%.1f min solver=%s",
                                    stage,
                                    String(profile),
                                    cand.id,
                                    elapsed_min,
                                    String(cand.values["solver_mode"])
                                )
                            )
                            next_heartbeat += TUNER_HEARTBEAT_SECONDS
                        end
                    end
                    wait(proc)
                    if !Base.success(proc)
                        run_ok = false
                        error_text = "nonzero_exit_code=$(proc.exitcode)"
                    end
                end
            catch err
                run_ok = false
                error_text = sprint(showerror, err)
            end
        end
    end

    summary_df = run_ok && isfile(paths.summary) ? CSV.read(paths.summary, DataFrame) : DataFrame()
    errors_df = isfile(paths.errors) ? CSV.read(paths.errors, DataFrame) : DataFrame()

    runtime_info = _parse_run_metrics(summary_df, wall_s)
    obj_rng = uncertainty_draw === nothing ? nothing : MersenneTwister(hash((cfg.seed, "telemetry_noise", cand.id, uncertainty_draw)))
    obj = _objective_from_summary(
        summary_df,
        cfg;
        run_failed=(!run_ok || nrow(summary_df) == 0),
        runtime_s=runtime_info.runtime_s,
        noise_rng=obj_rng
    )

    row = merge(
        _candidate_row_nt(cand),
        (
            stage=stage,
            profile=String(profile),
            uncertainty_draw=uncertainty_draw === nothing ? missing : uncertainty_draw,
            success=run_ok && nrow(summary_df) > 0,
            pass=obj.all_pass,
            objective=obj.objective,
            objective_base=obj.base_loss,
            penalty_fail=obj.fail_penalty,
            penalty_time=obj.runtime_penalty,
            failed_rows=obj.failed_rows,
            scenario_objective_breakdown=obj.per_scenario,
            scenario_count=runtime_info.scenario_count,
            event_count=runtime_info.event_count,
            odyssey_peri_rmse_km=_event_metric(summary_df, "odyssey", "peri", :rmse_km),
            odyssey_apo_rmse_km=_event_metric(summary_df, "odyssey", "apo", :rmse_km),
            vex_peri_rmse_km=_event_metric(summary_df, "vex", "peri", :rmse_km),
            vex_apo_rmse_km=_event_metric(summary_df, "vex", "apo", :rmse_km),
            earth_peri_rmse_km=_event_metric(summary_df, "earth_gmat", "peri", :rmse_km),
            earth_apo_rmse_km=_event_metric(summary_df, "earth_gmat", "apo", :rmse_km),
            simulation_runtime_s=runtime_info.runtime_s,
            wall_runtime_s=wall_s,
            error_message=run_ok ? "" : (isempty(error_text) ? "run_failed_or_summary_missing" : error_text),
            summary_path=paths.summary,
            errors_path=paths.errors,
            manifest_path=paths.manifest,
            log_path=paths.log,
            errors_rows=nrow(errors_df)
        )
    )
    return row
end

function _evaluate_tasks_batch(
    cfg::TunerConfig,
    stage::String,
    profile::Symbol,
    tasks::Vector{EvalTask},
    base_manifest::Dict{String, Any}
)::DataFrame
    n = length(tasks)
    n == 0 && return DataFrame()

    rows = Vector{NamedTuple}(undef, n)
    progress_lock = ReentrantLock()

    worker_count = min(cfg.parallel_candidates, n)
    if worker_count <= 1
        for idx in 1:n
            t = tasks[idx]
            _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage, idx, n, t.candidate.id))
            row = evaluate_candidate(
                cfg,
                stage,
                profile,
                t.candidate,
                base_manifest;
                uncertainty_draw=t.draw,
                progress_lock=progress_lock
            )
            rows[idx] = row
            _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d done success=%s objective=%.6f", stage, idx, n, t.candidate.id, string(row.success), Float64(row.objective)))
        end
        return DataFrame(rows)
    end

    jobs = Channel{Int}(n)
    for i in 1:n
        put!(jobs, i)
    end
    close(jobs)

    @sync for _ in 1:worker_count
        Threads.@spawn begin
            for idx in jobs
                t = tasks[idx]
                _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage, idx, n, t.candidate.id))
                row = evaluate_candidate(
                    cfg,
                    stage,
                    profile,
                    t.candidate,
                    base_manifest;
                    uncertainty_draw=t.draw,
                    progress_lock=progress_lock
                )
                rows[idx] = row
                _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d done success=%s objective=%.6f", stage, idx, n, t.candidate.id, string(row.success), Float64(row.objective)))
            end
        end
    end

    return DataFrame(rows)
end

function _candidate_from_row(row)::TuneCandidate
    values = Dict{String, Any}(
        "epoch_shift_s" => Int(row.epoch_shift_s),
        "ra_scale" => Float64(row.ra_scale),
        "rp_altitude_offset_m" => Float64(row.rp_altitude_offset_m),
        "i_offset_deg" => Float64(row.i_offset_deg),
        "aop_offset_deg" => Float64(row.aop_offset_deg),
        "raan_offset_deg" => Float64(row.raan_offset_deg),
        "ta_offset_deg" => Float64(row.ta_offset_deg),
        "bus_mass_scale" => Float64(row.bus_mass_scale),
        "prop_mass_scale" => Float64(row.prop_mass_scale),
        "panel_mass_scale" => Float64(row.panel_mass_scale),
        "bus_dims_scale" => Float64(row.bus_dims_scale),
        "panel_dims_scale" => Float64(row.panel_dims_scale),
        "panel_offset_scale" => Float64(row.panel_offset_scale),
        "srp_cr_scale" => Float64(row.srp_cr_scale),
        "srp_area_scale" => Float64(row.srp_area_scale),
        "gravity_variant" => String(row.gravity_variant),
        "dv_global_scale" => Float64(row.dv_global_scale),
        "dv_early_scale" => Float64(row.dv_early_scale),
        "dv_orbit7_bias_mps" => Float64(row.dv_orbit7_bias_mps),
        "solver_mode" => String(row.solver_mode_requested),
        "dt_max_orbit_s" => Float64(row.dt_max_orbit_requested_s),
        "dt_max_atm_s" => Float64(row.dt_max_atm_requested_s)
    )
    return TuneCandidate(id=Int(row.candidate_id), values=values, stage=String(row.stage))
end

function _sort_by_objective(df::DataFrame)::DataFrame
    good = df[df.success .== true, :]
    nrow(good) == 0 && return DataFrame()
    sort!(good, [:pass, :objective], rev=[true, false])
    return good
end

function _p95(vals::Vector{Float64})::Float64
    isempty(vals) && return Inf
    sorted = sort(vals)
    idx = max(1, ceil(Int, 0.95 * length(sorted)))
    return sorted[idx]
end

function _run_global_bo(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    base_manifest::Dict{String, Any}
)
    candidates = TuneCandidate[]
    scores = Float64[]
    all_rows = DataFrame()

    seen = Set{String}()
    next_id = 1

    init_candidates = _initial_design(cfg, params, cfg.n_init; start_id=next_id, stage="quick_global_init")
    next_id += length(init_candidates)
    append!(seen, (_candidate_signature(c) for c in init_candidates))

    init_tasks = [EvalTask(candidate=c) for c in init_candidates]
    init_df = _evaluate_tasks_batch(cfg, "quick_global_init", :quick, init_tasks, base_manifest)
    all_rows = vcat(all_rows, init_df; cols=:union)

    for row in eachrow(init_df)
        push!(candidates, _candidate_from_row(row))
        push!(scores, Float64(row.objective))
    end

    for k in 1:cfg.n_bo_iters
        model = _fit_surrogate(cfg, params, candidates, scores)
        model === nothing && break

        batch = _propose_batch(
            cfg,
            params,
            model,
            scores,
            seen,
            next_id,
            k,
            cfg.batch_size
        )
        isempty(batch) && break
        next_id += length(batch)

        stage = @sprintf("quick_global_bo_%03d", k)
        for i in eachindex(batch)
            batch[i] = TuneCandidate(id=batch[i].id, values=batch[i].values, stage=stage)
        end
        df = _evaluate_tasks_batch(cfg, stage, :quick, [EvalTask(candidate=c) for c in batch], base_manifest)
        all_rows = vcat(all_rows, df; cols=:union)

        for row in eachrow(df)
            push!(candidates, _candidate_from_row(row))
            push!(scores, Float64(row.objective))
        end
    end

    return all_rows, candidates, scores, next_id
end

function _continuous_params(params::Vector{ParameterSpec})::Vector{ParameterSpec}
    return [p for p in params if p.kind == :continuous]
end

function _run_local_refinement(
    cfg::TunerConfig,
    params::Vector{ParameterSpec},
    base_manifest::Dict{String, Any},
    quick_df::DataFrame,
    next_id::Int
)
    sorted = _sort_by_objective(quick_df)
    nrow(sorted) == 0 && return DataFrame(), next_id

    seed_count = min(cfg.local_seed_topk, nrow(sorted))
    seeds = [_candidate_from_row(sorted[i, :]) for i in 1:seed_count]
    cont = _continuous_params(params)
    isempty(cont) && return DataFrame(), next_id

    local_rows = DataFrame()
    stage_prefix = "quick_local"

    for (seed_idx, seed) in enumerate(seeds)
        center = seed
        center_score = Float64(sorted[seed_idx, :objective])
        step_scale = cfg.local_step_init

        for step in 1:cfg.local_steps
            neighbors = TuneCandidate[]
            for p in cont
                base_val = Float64(center.values[p.name])
                span = p.upper - p.lower
                delta = span * step_scale
                for sgn in (-1.0, 1.0)
                    vals = deepcopy(center.values)
                    vals[p.name] = clamp(base_val + sgn * delta, p.lower, p.upper)
                    cand = TuneCandidate(
                        id=next_id,
                        values=vals,
                        stage=@sprintf("%s_seed%02d_step%02d", stage_prefix, seed_idx, step)
                    )
                    push!(neighbors, cand)
                    next_id += 1
                end
            end

            isempty(neighbors) && break
            df = _evaluate_tasks_batch(cfg, neighbors[1].stage, :quick, [EvalTask(candidate=c) for c in neighbors], base_manifest)
            local_rows = vcat(local_rows, df; cols=:union)

            valid = df[df.success .== true, :]
            if nrow(valid) == 0
                step_scale *= cfg.local_step_shrink
                continue
            end
            sort!(valid, :objective)
            best_row = valid[1, :]
            best_score = Float64(best_row.objective)
            if best_score + cfg.local_min_improve < center_score
                center = _candidate_from_row(best_row)
                center_score = best_score
                step_scale *= cfg.local_step_expand
            else
                step_scale *= cfg.local_step_shrink
            end
        end
    end

    return local_rows, next_id
end

function _run_robustness_validation(
    cfg::TunerConfig,
    base_manifest::Dict{String, Any},
    finalists_df::DataFrame
)
    robust_samples = DataFrame()
    robust_rank_rows = NamedTuple[]

    for i in 1:nrow(finalists_df)
        cand = _candidate_from_row(finalists_df[i, :])
        tasks = [EvalTask(candidate=TuneCandidate(id=cand.id, values=deepcopy(cand.values), stage="robust_mc"), draw=j) for j in 1:cfg.robust_draws]
        sample_df = _evaluate_tasks_batch(cfg, "robust_mc", :quick, tasks, base_manifest)
        robust_samples = vcat(robust_samples, sample_df; cols=:union)

        success_mask = sample_df.success .== true
        sample_scores = Float64[sample_df.objective[j] for j in 1:nrow(sample_df) if success_mask[j]]
        fail_rate = nrow(sample_df) == 0 ? 1.0 : count(!, Bool.(sample_df.success)) / nrow(sample_df)

        if isempty(sample_scores)
            mean_j = Inf
            p95_j = Inf
            robust_j = Inf
        else
            mean_j = mean(sample_scores)
            p95_j = _p95(sample_scores)
            robust_j = mean_j + cfg.robust_alpha * p95_j + cfg.robust_beta * fail_rate
        end

        push!(robust_rank_rows, merge(
            _candidate_row_nt(cand),
            (
                stage="robust_rank",
                draws=cfg.robust_draws,
                mean_j=mean_j,
                p95_j=p95_j,
                fail_rate=fail_rate,
                j_robust=robust_j
            )
        ))
    end

    robust_rank_df = DataFrame(robust_rank_rows)
    if nrow(robust_rank_df) > 0
        sort!(robust_rank_df, :j_robust)
    end
    return robust_samples, robust_rank_df
end

function _write_report(
    path::String,
    cfg::TunerConfig,
    quick_df::DataFrame,
    robust_rank_df::DataFrame,
    full_df::DataFrame,
    best_manifest_path::String
)
    open(path, "w") do io
        println(io, "# Telemetry Hybrid Tuner Report")
        println(io)
        println(io, "- Generated (UTC): $(now(UTC))")
        println(io, "- Base manifest: `$(cfg.base_manifest)`")
        println(io, "- Outdir: `$(cfg.outdir)`")
        println(io, "- Global design: `n_init=$(cfg.n_init), n_bo_iters=$(cfg.n_bo_iters), batch=$(cfg.batch_size), acq=$(String(cfg.acquisition))`")
        println(io, "- Local refine: `topk=$(cfg.local_seed_topk), steps=$(cfg.local_steps)`")
        println(io, "- Robustness: `draws=$(cfg.robust_draws), alpha=$(cfg.robust_alpha), beta=$(cfg.robust_beta)`")
        println(io, "- Objective weights: `lambda_fail=$(cfg.lambda_fail), lambda_time=$(cfg.lambda_time), huber_delta=$(cfg.huber_delta)`")
        println(io)

        println(io, "## Objective")
        println(io)
        println(io, "`J(θ) = Σ_s w_s Σ_e [Huber(rmse_km/limit_rmse) + 0.5*Huber(max_abs_km/limit_abs)] + λfail*I(run_failed) + λtime*max(0, runtime/t_budget - 1)`")
        println(io)

        println(io, "## Stage Counts")
        println(io)
        println(io, "- Quick evaluations: `$(nrow(quick_df))`")
        println(io, "- Robustness samples: `$(nrow(robust_rank_df) == 0 ? 0 : nrow(robust_rank_df) * cfg.robust_draws)`")
        println(io, "- Full confirmations: `$(nrow(full_df))`")
        println(io)

        if nrow(robust_rank_df) > 0
            best = robust_rank_df[1, :]
            println(io, "## Best Robust Candidate")
            println(io)
            println(io, "- Candidate ID: `$(best.candidate_id)`")
            println(io, "- J_robust: `$(best.j_robust)`")
            println(io, "- Mean(J): `$(best.mean_j)`")
            println(io, "- P95(J): `$(best.p95_j)`")
            println(io, "- Fail rate: `$(best.fail_rate)`")
            println(io)
        end

        if nrow(full_df) > 0
            best_full = full_df[1, :]
            println(io, "## Full Profile Confirmation")
            println(io, "- Success: `$(best_full.success)`")
            println(io, "- Objective: `$(best_full.objective)`")
            println(io, "- Pass: `$(best_full.pass)`")
            println(io)
        end

        println(io, "## Artifacts")
        println(io)
        println(io, "- Quick results: `$(joinpath(cfg.outdir, "telemetry_hybrid_quick_results.csv"))`")
        println(io, "- Robust sample results: `$(joinpath(cfg.outdir, "telemetry_hybrid_robust_samples.csv"))`")
        println(io, "- Robust ranking: `$(joinpath(cfg.outdir, "telemetry_hybrid_robust_rank.csv"))`")
        println(io, "- Full results: `$(joinpath(cfg.outdir, "telemetry_hybrid_full_results.csv"))`")
        println(io, "- Combined results: `$(joinpath(cfg.outdir, "telemetry_hybrid_all_results.csv"))`")
        println(io, "- Best manifest: `$(best_manifest_path)`")
        println(io)

        println(io, "## Reproduce Best (Full)")
        println(io)
        println(io, "```bash")
        println(io, "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=.AGORA --startup-file=no src/analysis/verification/telemetry_verification.jl --profile=full --manifest=$(best_manifest_path) --enforce=$(cfg.enforce ? "1" : "0")")
        println(io, "```")
    end
end

function main_hybrid_tuner()
    cfg = parse_cli(copy(ARGS))
    mkpath(cfg.outdir)

    base_manifest = TOML.parsefile(cfg.base_manifest)
    params = _parameter_space(cfg)

    println("Telemetry hybrid tuner")
    println("outdir=$(cfg.outdir)")
    println("manifest=$(cfg.base_manifest)")
    println("global: n_init=$(cfg.n_init) n_bo_iters=$(cfg.n_bo_iters) batch=$(cfg.batch_size) acq=$(cfg.acquisition)")

    quick_global_df, _, _, next_id = _run_global_bo(cfg, params, base_manifest)

    println("running local refinement...")
    local_df, next_id = _run_local_refinement(cfg, params, base_manifest, quick_global_df, next_id)

    quick_df = vcat(quick_global_df, local_df; cols=:union)
    quick_csv = joinpath(cfg.outdir, "telemetry_hybrid_quick_results.csv")
    CSV.write(quick_csv, quick_df)

    valid_quick = _sort_by_objective(quick_df)
    nrow(valid_quick) > 0 || error("No successful quick evaluations produced by hybrid tuner.")

    finalist_count = min(cfg.finalists, nrow(valid_quick))
    finalists_df = valid_quick[1:finalist_count, :]

    println("running robustness Monte Carlo on $(finalist_count) finalists...")
    robust_samples_df, robust_rank_df = _run_robustness_validation(cfg, base_manifest, finalists_df)
    robust_samples_csv = joinpath(cfg.outdir, "telemetry_hybrid_robust_samples.csv")
    robust_rank_csv = joinpath(cfg.outdir, "telemetry_hybrid_robust_rank.csv")
    CSV.write(robust_samples_csv, robust_samples_df)
    CSV.write(robust_rank_csv, robust_rank_df)

    best_row = if nrow(robust_rank_df) > 0
        robust_rank_df[1, :]
    else
        finalists_df[1, :]
    end
    best_candidate = _candidate_from_row(best_row)

    full_df = DataFrame()
    if cfg.full_confirm
        println("confirming best robust candidate on full profile...")
        full_task = [EvalTask(candidate=TuneCandidate(id=best_candidate.id, values=deepcopy(best_candidate.values), stage="full_confirm"))]
        full_df = _evaluate_tasks_batch(cfg, "full_confirm", :full, full_task, base_manifest)
        sort!(full_df, :objective)
    end
    full_csv = joinpath(cfg.outdir, "telemetry_hybrid_full_results.csv")
    CSV.write(full_csv, full_df)

    best_manifest = _apply_candidate_to_manifest(base_manifest, best_candidate, cfg)
    best_manifest_path = joinpath(cfg.outdir, "telemetry_hybrid_best_manifest.toml")
    open(best_manifest_path, "w") do io
        TOML.print(io, best_manifest)
    end

    all_df = vcat(quick_df, robust_samples_df, full_df; cols=:union)
    all_csv = joinpath(cfg.outdir, "telemetry_hybrid_all_results.csv")
    CSV.write(all_csv, all_df)

    report_path = joinpath(cfg.outdir, "telemetry_hybrid_report.md")
    _write_report(report_path, cfg, quick_df, robust_rank_df, full_df, best_manifest_path)

    println()
    println("Hybrid tuning complete.")
    println("Quick results      : $quick_csv")
    println("Robust samples     : $robust_samples_csv")
    println("Robust rank        : $robust_rank_csv")
    println("Full results       : $full_csv")
    println("All results        : $all_csv")
    println("Best manifest      : $best_manifest_path")
    println("Report             : $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_hybrid_tuner()
end
