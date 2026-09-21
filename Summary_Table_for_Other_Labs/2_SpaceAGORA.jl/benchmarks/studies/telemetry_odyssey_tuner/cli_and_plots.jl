const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_BASE_MANIFEST = joinpath(REPO_ROOT, "test", "telemetry_benchmark_manifest.toml")
const DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "telemetry_tuning", "odyssey")
const TELEMETRY_STUDY_SCRIPT = joinpath(REPO_ROOT, "benchmarks", "studies", "telemetry_orbit_accuracy_study.jl")
const TUNER_HEARTBEAT_SECONDS = 120.0
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CSV
using DataFrames
using Dates
using Printf
using Random
using Statistics
using TOML

const PAPER_LIMITS = Dict(
    "peri" => (rmse_km=1.715, max_abs_km=8.4),
    "apo" => (rmse_km=653.057, max_abs_km=1109.86)
)
const _PLOTS_READY = Ref(false)
const _PLOTS_AVAILABLE = Ref(false)

const EPOCH_SHIFT_VALUES_S = collect(-180:60:180)
const RA_SCALE_VALUES = [0.999, 1.0, 1.001]
const RP_ALTITUDE_OFFSET_VALUES_M = [-3_000.0, 0.0, 3_000.0]
const I_OFFSET_VALUES_DEG = [-0.03, 0.0, 0.03]
const AOP_OFFSET_VALUES_DEG = [-0.2, 0.0, 0.2]
const RAAN_OFFSET_VALUES_DEG = [-0.2, 0.0, 0.2]
const TA_OFFSET_VALUES_DEG = [-0.5, 0.0, 0.5]
const BUS_MASS_SCALE_VALUES = [0.99, 1.0, 1.01]
const PROP_MASS_SCALE_VALUES = [0.95, 1.0, 1.05]
const PANEL_MASS_SCALE_VALUES = [0.95, 1.0, 1.05]
const BUS_DIMS_SCALE_VALUES = [0.995, 1.0, 1.005]
const PANEL_DIMS_SCALE_VALUES = [0.995, 1.0, 1.005]
const PANEL_OFFSET_SCALE_VALUES = [0.99, 1.0, 1.01]
const SRP_CR_SCALE_VALUES = [0.95, 1.0, 1.05]
const SRP_AREA_SCALE_VALUES = [0.95, 1.0, 1.05]
const GRAVITY_VARIANTS = ["harm_deg2", "harm_deg4", "harm_deg8", "harm_deg20", "j2_only"]
const DV_GLOBAL_SCALE_VALUES = [0.8, 1.0, 1.2]
const DV_EARLY_SCALE_VALUES = [0.8, 1.0, 1.2]
const DV_ORBIT7_BIAS_VALUES_MPS = [-0.3, 0.0, 0.3]

Base.@kwdef struct TunerConfig
    base_manifest::String = DEFAULT_BASE_MANIFEST
    outdir::String = DEFAULT_OUTDIR
    scenario_name::String = "odyssey"
    quick_candidates::Int = 24
    full_topk::Int = 4
    parallel_candidates::Int = 1
    seed::Int = 42
    quick_only::Bool = false
    enforce::Bool = false
    calibration_enabled::Bool = false
    solver_modes::Vector{String} = String["auto_stiff"]
    dt_max_orbit_values::Vector{Float64} = [30.0, 60.0]
    dt_max_atm_values::Vector{Float64} = [0.1, 0.2]
    maxiters_quick::Union{Nothing, Int} = nothing
    maxiters_full::Union{Nothing, Int} = nothing
end

Base.@kwdef struct TuneCandidate
    id::Int
    epoch_shift_s::Int
    ra_scale::Float64
    rp_altitude_offset_m::Float64
    i_offset_deg::Float64
    aop_offset_deg::Float64
    raan_offset_deg::Float64
    ta_offset_deg::Float64
    bus_mass_scale::Float64
    prop_mass_scale::Float64
    panel_mass_scale::Float64
    bus_dims_scale::Float64
    panel_dims_scale::Float64
    panel_offset_scale::Float64
    srp_cr_scale::Float64
    srp_area_scale::Float64
    gravity_variant::String
    dv_global_scale::Float64
    dv_early_scale::Float64
    dv_orbit7_bias_mps::Float64
    solver_mode::String
    dt_max_orbit_s::Float64
    dt_max_atm_s::Float64
end

@inline function _parse_bool(raw::AbstractString)::Bool
    token = lowercase(strip(raw))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean token '$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _parse_int(raw::AbstractString, label::String)::Int
    value = try
        parse(Int, raw)
    catch
        throw(ArgumentError("Expected integer for $label, got '$raw'"))
    end
    return value
end

@inline function _parse_float(raw::AbstractString, label::String)::Float64
    value = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("Expected float for $label, got '$raw'"))
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
    isempty(vals) && throw(ArgumentError("Expected non-empty comma-separated float list for $label"))
    return vals
end

@inline function _parse_string_list(raw::AbstractString, label::String)::Vector{String}
    vals = String[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(vals, t)
    end
    isempty(vals) && throw(ArgumentError("Expected non-empty comma-separated list for $label"))
    return vals
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
    # Odyssey candidates are single-spacecraft and usually run in outer process-level
    # parallelism, so keep inner callback/effector threading disabled by default.
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

function _ensure_plots_loaded()::Bool
    if _PLOTS_READY[]
        return _PLOTS_AVAILABLE[]
    end
    _PLOTS_READY[] = true
    try
        @eval import Plots
        _PLOTS_AVAILABLE[] = true
    catch err
        @warn "Plots.jl is unavailable; tuner plots will be skipped." exception=(err, catch_backtrace())
        _PLOTS_AVAILABLE[] = false
    end
    return _PLOTS_AVAILABLE[]
end

function parse_tuner_cli(args::Vector{String})::TunerConfig
    base_manifest = DEFAULT_BASE_MANIFEST
    outdir = DEFAULT_OUTDIR
    scenario_name = "odyssey"
    quick_candidates = 24
    full_topk = 4
    parallel_candidates = 1
    seed = 42
    quick_only = false
    enforce = false
    calibration_enabled = false
    solver_modes = String["auto_stiff"]
    dt_max_orbit_values = [30.0, 60.0]
    dt_max_atm_values = [0.1, 0.2]
    maxiters_quick = nothing
    maxiters_full = nothing

    for arg in args
        if startswith(arg, "--base-manifest=")
            base_manifest = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--outdir=")
            outdir = abspath(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--scenario=")
            scenario_name = strip(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--quick-candidates=")
            quick_candidates = _parse_int(split(arg, "=", limit=2)[2], "quick-candidates")
        elseif startswith(arg, "--topk-full=")
            full_topk = _parse_int(split(arg, "=", limit=2)[2], "topk-full")
        elseif startswith(arg, "--parallel-candidates=")
            parallel_candidates = _parse_int(split(arg, "=", limit=2)[2], "parallel-candidates")
        elseif startswith(arg, "--seed=")
            seed = _parse_int(split(arg, "=", limit=2)[2], "seed")
        elseif startswith(arg, "--quick-only=")
            quick_only = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--enforce=")
            enforce = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--calibration=")
            calibration_enabled = _parse_bool(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--solver-modes=")
            solver_modes = _parse_string_list(split(arg, "=", limit=2)[2], "solver-modes")
        elseif startswith(arg, "--dt-orbit-values=")
            dt_max_orbit_values = _parse_float_list(split(arg, "=", limit=2)[2], "dt-orbit-values")
        elseif startswith(arg, "--dt-atm-values=")
            dt_max_atm_values = _parse_float_list(split(arg, "=", limit=2)[2], "dt-atm-values")
        elseif startswith(arg, "--maxiters-quick=")
            maxiters_quick = _parse_int(split(arg, "=", limit=2)[2], "maxiters-quick")
        elseif startswith(arg, "--maxiters-full=")
            maxiters_full = _parse_int(split(arg, "=", limit=2)[2], "maxiters-full")
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: --base-manifest=..., --outdir=..., --scenario=..., " *
                "--quick-candidates=N, --topk-full=K, --parallel-candidates=N, --seed=N, --quick-only=0|1, --enforce=0|1, " *
                "--calibration=0|1, --solver-modes=a,b, --dt-orbit-values=a,b, --dt-atm-values=a,b, " *
                "--maxiters-quick=N, --maxiters-full=N"
            ))
        end
    end
    quick_candidates > 0 || throw(ArgumentError("quick-candidates must be > 0"))
    full_topk > 0 || throw(ArgumentError("topk-full must be > 0"))
    parallel_candidates > 0 || throw(ArgumentError("parallel-candidates must be > 0"))
    all(v -> v > 0.0, dt_max_orbit_values) || throw(ArgumentError("dt-orbit-values must be > 0"))
    all(v -> v > 0.0, dt_max_atm_values) || throw(ArgumentError("dt-atm-values must be > 0"))
    return TunerConfig(
        base_manifest=base_manifest,
        outdir=outdir,
        scenario_name=scenario_name,
        quick_candidates=quick_candidates,
        full_topk=full_topk,
        parallel_candidates=parallel_candidates,
        seed=seed,
        quick_only=quick_only,
        enforce=enforce,
        calibration_enabled=calibration_enabled,
        solver_modes=solver_modes,
        dt_max_orbit_values=dt_max_orbit_values,
        dt_max_atm_values=dt_max_atm_values,
        maxiters_quick=maxiters_quick,
        maxiters_full=maxiters_full
    )
end

@inline function _wrap_0_360(deg::Float64)::Float64
    out = mod(deg, 360.0)
    out < 0.0 && return out + 360.0
    return out
end

@inline _clamp_positive(x::Float64) = max(1e-9, x)
@inline _to_f64(v) = v isa Missing ? NaN : Float64(v)
@inline function _plots_module()
    return Base.invokelatest(() -> getproperty(Main, :Plots))
end
@inline _plots_fn(sym::Symbol) = getproperty(_plots_module(), sym)
@inline function _plots_call(sym::Symbol, args...; kwargs...)
    return Base.invokelatest(_plots_fn(sym), args...; kwargs...)
end

function _plot_scatter_candidate_score(df::DataFrame, stage::String, path::String)::Bool
    nrow(df) == 0 && return false
    mask = [isfinite(_to_f64(df.score[i])) for i in 1:nrow(df)]
    any(mask) || return false
    ids = Int[df.candidate_id[i] for i in 1:nrow(df) if mask[i]]
    scores = Float64[_to_f64(df.score[i]) for i in 1:nrow(df) if mask[i]]
    passed = Bool[df.pass[i] for i in 1:nrow(df) if mask[i]]
    colors = [p ? :forestgreen : :firebrick for p in passed]
    marks = [p ? :circle : :x for p in passed]
    p = _plots_call(
        :scatter,
        ids,
        scores;
        xlabel="candidate_id",
        ylabel="score (lower is better)",
        title="$(uppercasefirst(stage)) Candidate Score",
        color=colors,
        markershape=marks,
        markerstrokewidth=0,
        legend=false,
        size=(1000, 650)
    )
    _plots_call(:savefig, p, path)
    return true
end

function _plot_tradeoff_peri_apo(df::DataFrame, stage::String, path::String)::Bool
    nrow(df) == 0 && return false
    mask = [
        isfinite(_to_f64(df.peri_rmse_km[i])) &&
        isfinite(_to_f64(df.apo_rmse_km[i]))
        for i in 1:nrow(df)
    ]
    any(mask) || return false
    peri = Float64[_to_f64(df.peri_rmse_km[i]) for i in 1:nrow(df) if mask[i]]
    apo = Float64[_to_f64(df.apo_rmse_km[i]) for i in 1:nrow(df) if mask[i]]
    passed = Bool[df.pass[i] for i in 1:nrow(df) if mask[i]]
    labels = [string(df.candidate_id[i]) for i in 1:nrow(df) if mask[i]]
    colors = [p ? :forestgreen : :firebrick for p in passed]

    p = _plots_call(
        :scatter,
        apo,
        peri;
        xlabel="apo RMSE [km]",
        ylabel="peri RMSE [km]",
        title="$(uppercasefirst(stage)) RMSE Tradeoff",
        color=colors,
        markerstrokewidth=0,
        legend=false,
        size=(1000, 650)
    )
    for i in eachindex(apo)
        txt = _plots_call(:text, labels[i], 7, :gray35)
        _plots_call(:annotate!, p, apo[i], peri[i], txt)
    end
    _plots_call(:vline!, p, [PAPER_LIMITS["apo"].rmse_km]; color=:black, linestyle=:dash, linewidth=1.0)
    _plots_call(:hline!, p, [PAPER_LIMITS["peri"].rmse_km]; color=:black, linestyle=:dash, linewidth=1.0)
    _plots_call(:savefig, p, path)
    return true
end

function _plot_topk_bar(df::DataFrame, stage::String, path::String; k::Int=10)::Bool
    nrow(df) == 0 && return false
    work = df[[isfinite(_to_f64(v)) for v in df.score], :]
    nrow(work) == 0 && return false
    sort!(work, :score)
    n = min(k, nrow(work))
    sub = work[1:n, :]
    labels = [string(sub.candidate_id[i]) for i in 1:n]
    vals = Float64[_to_f64(sub.score[i]) for i in 1:n]
    colors = [Bool(sub.pass[i]) ? :forestgreen : :firebrick for i in 1:n]
    p = _plots_call(
        :bar,
        labels,
        vals;
        xlabel="candidate_id",
        ylabel="score",
        title="$(uppercasefirst(stage)) Top $(n) Scores",
        color=colors,
        legend=false,
        size=(1100, 650)
    )
    _plots_call(:savefig, p, path)
    return true
end

function _generate_tuning_plots(outdir::String, quick_df::DataFrame, full_df::DataFrame)::Vector{String}
    _ensure_plots_loaded() || return String[]
    plots_dir = joinpath(outdir, "plots")
    mkpath(plots_dir)
    paths = String[]

    quick_score = joinpath(plots_dir, "quick_score_vs_candidate.png")
    if _plot_scatter_candidate_score(quick_df, "quick", quick_score)
        push!(paths, quick_score)
    end
    quick_tradeoff = joinpath(plots_dir, "quick_tradeoff_peri_vs_apo_rmse.png")
    if _plot_tradeoff_peri_apo(quick_df, "quick", quick_tradeoff)
        push!(paths, quick_tradeoff)
    end
    quick_top = joinpath(plots_dir, "quick_top_scores.png")
    if _plot_topk_bar(quick_df, "quick", quick_top; k=10)
        push!(paths, quick_top)
    end

    if nrow(full_df) > 0
        full_score = joinpath(plots_dir, "full_score_vs_candidate.png")
        if _plot_scatter_candidate_score(full_df, "full", full_score)
            push!(paths, full_score)
        end
        full_tradeoff = joinpath(plots_dir, "full_tradeoff_peri_vs_apo_rmse.png")
        if _plot_tradeoff_peri_apo(full_df, "full", full_tradeoff)
            push!(paths, full_tradeoff)
        end
        full_top = joinpath(plots_dir, "full_top_scores.png")
        if _plot_topk_bar(full_df, "full", full_top; k=10)
            push!(paths, full_top)
        end
    end

    return paths
end

@inline function _candidate_signature(c::TuneCandidate)::String
    return join((
        string(c.epoch_shift_s),
        @sprintf("%.6f", c.ra_scale),
        @sprintf("%.3f", c.rp_altitude_offset_m),
        @sprintf("%.5f", c.i_offset_deg),
        @sprintf("%.5f", c.aop_offset_deg),
        @sprintf("%.5f", c.raan_offset_deg),
        @sprintf("%.5f", c.ta_offset_deg),
        @sprintf("%.6f", c.bus_mass_scale),
        @sprintf("%.6f", c.prop_mass_scale),
        @sprintf("%.6f", c.panel_mass_scale),
        @sprintf("%.6f", c.bus_dims_scale),
        @sprintf("%.6f", c.panel_dims_scale),
        @sprintf("%.6f", c.panel_offset_scale),
        @sprintf("%.6f", c.srp_cr_scale),
        @sprintf("%.6f", c.srp_area_scale),
        c.gravity_variant,
        @sprintf("%.6f", c.dv_global_scale),
        @sprintf("%.6f", c.dv_early_scale),
        @sprintf("%.6f", c.dv_orbit7_bias_mps),
        c.solver_mode,
        @sprintf("%.3f", c.dt_max_orbit_s),
        @sprintf("%.3f", c.dt_max_atm_s)
    ), "|")
end
