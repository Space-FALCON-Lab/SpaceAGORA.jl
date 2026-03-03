const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_BASE_MANIFEST = joinpath(@__DIR__, "telemetry_benchmark_manifest.toml")
const DEFAULT_OUTDIR = joinpath(REPO_ROOT, "output", "telemetry_tuning", "odyssey")
const TELEMETRY_STUDY_SCRIPT = joinpath(REPO_ROOT, "scripts", "verify_telemetry.jl")
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

function _candidate_nt(c::TuneCandidate)
    return (
        candidate_id=c.id,
        epoch_shift_s=c.epoch_shift_s,
        ra_scale=c.ra_scale,
        rp_altitude_offset_m=c.rp_altitude_offset_m,
        i_offset_deg=c.i_offset_deg,
        aop_offset_deg=c.aop_offset_deg,
        raan_offset_deg=c.raan_offset_deg,
        ta_offset_deg=c.ta_offset_deg,
        bus_mass_scale=c.bus_mass_scale,
        prop_mass_scale=c.prop_mass_scale,
        panel_mass_scale=c.panel_mass_scale,
        bus_dims_scale=c.bus_dims_scale,
        panel_dims_scale=c.panel_dims_scale,
        panel_offset_scale=c.panel_offset_scale,
        srp_cr_scale=c.srp_cr_scale,
        srp_area_scale=c.srp_area_scale,
        gravity_variant=c.gravity_variant,
        dv_global_scale=c.dv_global_scale,
        dv_early_scale=c.dv_early_scale,
        dv_orbit7_bias_mps=c.dv_orbit7_bias_mps,
        solver_mode_requested=c.solver_mode,
        dt_max_orbit_requested_s=c.dt_max_orbit_s,
        dt_max_atm_requested_s=c.dt_max_atm_s
    )
end

@inline _pick(rng::AbstractRNG, values::Vector) = values[rand(rng, 1:length(values))]

function sample_candidates(cfg::TunerConfig)::Vector{TuneCandidate}
    rng = MersenneTwister(cfg.seed)
    candidates = TuneCandidate[]
    seen = Set{String}()
    while length(candidates) < cfg.quick_candidates
        cand = TuneCandidate(
            id=length(candidates) + 1,
            epoch_shift_s=_pick(rng, EPOCH_SHIFT_VALUES_S),
            ra_scale=_pick(rng, RA_SCALE_VALUES),
            rp_altitude_offset_m=_pick(rng, RP_ALTITUDE_OFFSET_VALUES_M),
            i_offset_deg=_pick(rng, I_OFFSET_VALUES_DEG),
            aop_offset_deg=_pick(rng, AOP_OFFSET_VALUES_DEG),
            raan_offset_deg=_pick(rng, RAAN_OFFSET_VALUES_DEG),
            ta_offset_deg=_pick(rng, TA_OFFSET_VALUES_DEG),
            bus_mass_scale=_pick(rng, BUS_MASS_SCALE_VALUES),
            prop_mass_scale=_pick(rng, PROP_MASS_SCALE_VALUES),
            panel_mass_scale=_pick(rng, PANEL_MASS_SCALE_VALUES),
            bus_dims_scale=_pick(rng, BUS_DIMS_SCALE_VALUES),
            panel_dims_scale=_pick(rng, PANEL_DIMS_SCALE_VALUES),
            panel_offset_scale=_pick(rng, PANEL_OFFSET_SCALE_VALUES),
            srp_cr_scale=_pick(rng, SRP_CR_SCALE_VALUES),
            srp_area_scale=_pick(rng, SRP_AREA_SCALE_VALUES),
            gravity_variant=_pick(rng, GRAVITY_VARIANTS),
            dv_global_scale=_pick(rng, DV_GLOBAL_SCALE_VALUES),
            dv_early_scale=_pick(rng, DV_EARLY_SCALE_VALUES),
            dv_orbit7_bias_mps=_pick(rng, DV_ORBIT7_BIAS_VALUES_MPS),
            solver_mode=_pick(rng, cfg.solver_modes),
            dt_max_orbit_s=_pick(rng, cfg.dt_max_orbit_values),
            dt_max_atm_s=_pick(rng, cfg.dt_max_atm_values)
        )
        sig = _candidate_signature(cand)
        if sig in seen
            continue
        end
        push!(seen, sig)
        push!(candidates, cand)
    end
    return candidates
end

function _find_base_scenario(base_manifest::Dict{String, Any}, name::String)::Dict{String, Any}
    haskey(base_manifest, "scenarios") || throw(ArgumentError("Manifest is missing top-level 'scenarios' array"))
    scenarios = base_manifest["scenarios"]
    scenarios isa AbstractVector || throw(ArgumentError("'scenarios' must be an array in manifest"))
    for sc in scenarios
        if sc isa AbstractDict && haskey(sc, "name") && String(sc["name"]) == name
            return deepcopy(Dict{String, Any}(sc))
        end
    end
    throw(ArgumentError("Scenario '$name' not found in manifest $(DEFAULT_BASE_MANIFEST)"))
end

function _apply_candidate(base::Dict{String, Any}, cand::TuneCandidate; calibration_enabled::Bool)::Dict{String, Any}
    sc = deepcopy(base)

    init = Dict{String, Any}(sc["initial_time"])
    sec_raw = Float64(init["second"])
    sec_int = floor(Int, sec_raw)
    sec_frac = sec_raw - sec_int
    sec_ms = round(Int, sec_frac * 1000)
    dt0 = DateTime(Int(init["year"]), Int(init["month"]), Int(init["day"]), Int(init["hour"]), Int(init["minute"]), sec_int, sec_ms)
    dt_shifted = dt0 + Second(cand.epoch_shift_s)
    init["year"] = year(dt_shifted)
    init["month"] = month(dt_shifted)
    init["day"] = day(dt_shifted)
    init["hour"] = hour(dt_shifted)
    init["minute"] = minute(dt_shifted)
    init["second"] = second(dt_shifted) + millisecond(dt_shifted) / 1000
    sc["initial_time"] = init

    sc["ra_m"] = _clamp_positive(Float64(sc["ra_m"]) * cand.ra_scale)
    sc["rp_altitude_m"] = _clamp_positive(Float64(sc["rp_altitude_m"]) + cand.rp_altitude_offset_m)
    sc["i_deg"] = clamp(Float64(sc["i_deg"]) + cand.i_offset_deg, 0.0, 180.0)
    sc["aop_deg"] = _wrap_0_360(Float64(sc["aop_deg"]) + cand.aop_offset_deg)
    sc["raan_deg"] = _wrap_0_360(Float64(sc["raan_deg"]) + cand.raan_offset_deg)
    sc["ta_deg"] = _wrap_0_360(Float64(sc["ta_deg"]) + cand.ta_offset_deg)

    if haskey(sc, "spacecraft")
        spacecraft = Dict{String, Any}(sc["spacecraft"])
        spacecraft["bus_mass_kg"] = _clamp_positive(Float64(spacecraft["bus_mass_kg"]) * cand.bus_mass_scale)
        spacecraft["prop_mass_kg"] = _clamp_positive(Float64(spacecraft["prop_mass_kg"]) * cand.prop_mass_scale)
        spacecraft["panel_mass_each_kg"] = _clamp_positive(Float64(spacecraft["panel_mass_each_kg"]) * cand.panel_mass_scale)
        spacecraft["panel_offset_y_m"] = _clamp_positive(Float64(spacecraft["panel_offset_y_m"]) * cand.panel_offset_scale)
        spacecraft["bus_dims_m"] = [Float64(v) * cand.bus_dims_scale for v in spacecraft["bus_dims_m"]]
        spacecraft["panel_dims_m"] = [Float64(v) * cand.panel_dims_scale for v in spacecraft["panel_dims_m"]]
        sc["spacecraft"] = spacecraft
    end

    if haskey(sc, "srp_cr")
        sc["srp_cr"] = _clamp_positive(Float64(sc["srp_cr"]) * cand.srp_cr_scale)
    end
    if haskey(sc, "srp_area_m2")
        sc["srp_area_m2"] = _clamp_positive(Float64(sc["srp_area_m2"]) * cand.srp_area_scale)
    end

    if cand.gravity_variant == "j2_only"
        sc["gravity_model"] = "inverse_squared_j2"
        sc["gravity_harmonics_degree"] = 0
        sc["gravity_harmonics_order"] = 0
    else
        degree = if cand.gravity_variant == "harm_deg2"
            2
        elseif cand.gravity_variant == "harm_deg4"
            4
        elseif cand.gravity_variant == "harm_deg8"
            8
        else
            20
        end
        sc["gravity_model"] = "inverse_squared"
        sc["gravity_harmonics_degree"] = degree
        sc["gravity_harmonics_order"] = degree
    end

    if haskey(sc, "maneuvers")
        maneuvers = Dict{String, Any}(sc["maneuvers"])
        orbits = [Int(v) for v in maneuvers["orbit_numbers"]]
        delta_v = [Float64(v) for v in maneuvers["delta_v_mps"]]
        tuned = similar(delta_v)
        for i in eachindex(delta_v)
            scale = cand.dv_global_scale
            if orbits[i] <= 50
                scale *= cand.dv_early_scale
            end
            dv = delta_v[i] * scale
            if orbits[i] == 7
                dv += cand.dv_orbit7_bias_mps
            end
            tuned[i] = dv
        end
        maneuvers["delta_v_mps"] = tuned
        sc["maneuvers"] = maneuvers
    end

    if haskey(sc, "calibration")
        cal = Dict{String, Any}(sc["calibration"])
        cal["enabled"] = calibration_enabled
        sc["calibration"] = cal
    end
    return sc
end

@inline function _paper_ratio(event::String, rmse_km::Float64, max_abs_km::Float64)
    lim = PAPER_LIMITS[event]
    return rmse_km / lim.rmse_km, max_abs_km / lim.max_abs_km
end

@inline function _score_candidate(
    peri_rmse::Float64,
    apo_rmse::Float64,
    peri_abs::Float64,
    apo_abs::Float64
)::Float64
    peri_rmse_ratio, peri_abs_ratio = _paper_ratio("peri", peri_rmse, peri_abs)
    apo_rmse_ratio, apo_abs_ratio = _paper_ratio("apo", apo_rmse, apo_abs)
    ratios = [peri_rmse_ratio, peri_abs_ratio, apo_rmse_ratio, apo_abs_ratio]
    return mean(ratios) + 0.35 * maximum(ratios)
end

@inline function _passes_paper_limits(
    peri_rmse::Float64,
    apo_rmse::Float64,
    peri_abs::Float64,
    apo_abs::Float64
)::Bool
    peri_ok = peri_rmse <= PAPER_LIMITS["peri"].rmse_km && peri_abs <= PAPER_LIMITS["peri"].max_abs_km
    apo_ok = apo_rmse <= PAPER_LIMITS["apo"].rmse_km && apo_abs <= PAPER_LIMITS["apo"].max_abs_km
    return peri_ok && apo_ok
end

function _event_row(df::DataFrame, event::String)
    idx = findfirst((df.event .== event))
    idx === nothing && throw(ArgumentError("Missing event '$event' in summary dataframe"))
    return df[idx, :]
end

function _candidate_paths(outdir::String, stage::Symbol, id::Int)
    tag = @sprintf("cand_%04d_%s", id, String(stage))
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

function evaluate_candidate(
    cfg::TunerConfig,
    stage::Symbol,
    cand::TuneCandidate,
    tuned_scenario::Dict{String, Any};
    progress_lock::Union{Nothing, ReentrantLock}=nothing
)
    paths = _candidate_paths(cfg.outdir, stage, cand.id)
    manifest_doc = Dict{String, Any}("version" => 1, "scenarios" => Any[tuned_scenario])
    open(paths.manifest, "w") do io
        TOML.print(io, manifest_doc)
    end

    profile = String(stage)
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(REPO_ROOT, ".AGORA")) --startup-file=no $(TELEMETRY_STUDY_SCRIPT) --profile=$(profile) --manifest=$(paths.manifest) --enforce=$(cfg.enforce ? "1" : "0") --out-summary=$(paths.summary) --out-errors=$(paths.errors)`

    env_pairs = Pair{String, String}[
        "SPACEAGORA_TELEMETRY_SOLVER_MODE" => cand.solver_mode,
        "SPACEAGORA_TELEMETRY_DT_MAX_ORBIT" => string(cand.dt_max_orbit_s),
        "SPACEAGORA_TELEMETRY_DT_MAX_ATM" => string(cand.dt_max_atm_s)
    ]
    append!(env_pairs, _candidate_runtime_policy_env_pairs())
    if stage == :quick && cfg.maxiters_quick !== nothing
        push!(env_pairs, "SPACEAGORA_TELEMETRY_SOLVER_MAXITERS" => string(cfg.maxiters_quick))
    elseif stage == :full && cfg.maxiters_full !== nothing
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
                                    "[heartbeat %s] candidate=%d elapsed=%.1f min solver=%s dt_orbit=%.1f dt_atm=%.3f",
                                    profile,
                                    cand.id,
                                    elapsed_min,
                                    cand.solver_mode,
                                    cand.dt_max_orbit_s,
                                    cand.dt_max_atm_s
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

    if !run_ok || !isfile(paths.summary)
        return merge(
            _candidate_nt(cand),
            (
                stage=profile,
                success=false,
                pass=false,
                score=Inf,
                error_message=isempty(error_text) ? "run_failed_or_summary_missing" : error_text,
                peri_rmse_km=missing,
                apo_rmse_km=missing,
                peri_max_abs_km=missing,
                apo_max_abs_km=missing,
                peri_nmae=missing,
                apo_nmae=missing,
                solver_mode_reported=missing,
                solver_sequence=missing,
                solver_fallback_count=missing,
                solver_retcode=missing,
                coverage_peri=missing,
                coverage_apo=missing,
                simulation_runtime_s=missing,
                wall_runtime_s=wall_s,
                summary_path=paths.summary,
                errors_path=paths.errors,
                manifest_path=paths.manifest,
                log_path=paths.log
            )
        )
    end

    summary_df = CSV.read(paths.summary, DataFrame)
    scenario_df = summary_df[summary_df.scenario .== cfg.scenario_name, :]
    if nrow(scenario_df) == 0
        return merge(
            _candidate_nt(cand),
            (
                stage=profile,
                success=false,
                pass=false,
                score=Inf,
                error_message="summary_missing_scenario",
                peri_rmse_km=missing,
                apo_rmse_km=missing,
                peri_max_abs_km=missing,
                apo_max_abs_km=missing,
                peri_nmae=missing,
                apo_nmae=missing,
                solver_mode_reported=missing,
                solver_sequence=missing,
                solver_fallback_count=missing,
                solver_retcode=missing,
                coverage_peri=missing,
                coverage_apo=missing,
                simulation_runtime_s=missing,
                wall_runtime_s=wall_s,
                summary_path=paths.summary,
                errors_path=paths.errors,
                manifest_path=paths.manifest,
                log_path=paths.log
            )
        )
    end

    peri = _event_row(scenario_df, "peri")
    apo = _event_row(scenario_df, "apo")
    peri_rmse = Float64(peri.rmse_km)
    apo_rmse = Float64(apo.rmse_km)
    peri_abs = Float64(peri.max_abs_km)
    apo_abs = Float64(apo.max_abs_km)
    score = _score_candidate(peri_rmse, apo_rmse, peri_abs, apo_abs)
    pass = _passes_paper_limits(peri_rmse, apo_rmse, peri_abs, apo_abs)

    return merge(
        _candidate_nt(cand),
        (
            stage=profile,
            success=true,
            pass=pass,
            score=score,
            error_message="",
            peri_rmse_km=peri_rmse,
            apo_rmse_km=apo_rmse,
            peri_max_abs_km=peri_abs,
            apo_max_abs_km=apo_abs,
            peri_nmae=Float64(peri.nmae),
            apo_nmae=Float64(apo.nmae),
            solver_mode_reported=String(peri.solver_mode),
            solver_sequence=String(peri.solver_sequence),
            solver_fallback_count=Int(peri.solver_fallback_count),
            solver_retcode=String(peri.solver_retcode),
            coverage_peri=Float64(peri.coverage),
            coverage_apo=Float64(apo.coverage),
            simulation_runtime_s=Float64(peri.simulation_runtime_s),
            wall_runtime_s=wall_s,
            summary_path=paths.summary,
            errors_path=paths.errors,
            manifest_path=paths.manifest,
            log_path=paths.log
        )
    )
end

function _evaluate_candidates_batch(
    cfg::TunerConfig,
    stage::Symbol,
    candidates::Vector{TuneCandidate},
    base_scenario::Dict{String, Any}
)::DataFrame
    n = length(candidates)
    n == 0 && return DataFrame()
    stage_name = String(stage)
    rows = Vector{NamedTuple}(undef, n)
    progress_lock = ReentrantLock()

    worker_count = min(cfg.parallel_candidates, n)
    if worker_count <= 1
        for idx in 1:n
            cand = candidates[idx]
            _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage_name, idx, n, cand.id))
            tuned = _apply_candidate(base_scenario, cand; calibration_enabled=cfg.calibration_enabled)
            row = evaluate_candidate(cfg, stage, cand, tuned; progress_lock=progress_lock)
            rows[idx] = row
            _progress_print(
                progress_lock,
                @sprintf(
                    "[%s %d/%d] candidate=%d done success=%s pass=%s score=%s wall_s=%.1f",
                    stage_name,
                    idx,
                    n,
                    cand.id,
                    string(row.success),
                    string(row.pass),
                    string(row.score),
                    Float64(row.wall_runtime_s)
                )
            )
        end
        return DataFrame(rows)
    end

    jobs = Channel{Int}(n)
    for idx in 1:n
        put!(jobs, idx)
    end
    close(jobs)

    @sync for _ in 1:worker_count
        Threads.@spawn begin
            for idx in jobs
                cand = candidates[idx]
                _progress_print(progress_lock, @sprintf("[%s %d/%d] candidate=%d start", stage_name, idx, n, cand.id))
                tuned = _apply_candidate(base_scenario, cand; calibration_enabled=cfg.calibration_enabled)
                row = evaluate_candidate(cfg, stage, cand, tuned; progress_lock=progress_lock)
                rows[idx] = row
                _progress_print(
                    progress_lock,
                    @sprintf(
                        "[%s %d/%d] candidate=%d done success=%s pass=%s score=%s wall_s=%.1f",
                        stage_name,
                        idx,
                        n,
                        cand.id,
                        string(row.success),
                        string(row.pass),
                        string(row.score),
                        Float64(row.wall_runtime_s)
                    )
                )
            end
        end
    end
    return DataFrame(rows)
end

function _candidate_from_row(row)::TuneCandidate
    return TuneCandidate(
        id=Int(row.candidate_id),
        epoch_shift_s=Int(row.epoch_shift_s),
        ra_scale=Float64(row.ra_scale),
        rp_altitude_offset_m=Float64(row.rp_altitude_offset_m),
        i_offset_deg=Float64(row.i_offset_deg),
        aop_offset_deg=Float64(row.aop_offset_deg),
        raan_offset_deg=Float64(row.raan_offset_deg),
        ta_offset_deg=Float64(row.ta_offset_deg),
        bus_mass_scale=Float64(row.bus_mass_scale),
        prop_mass_scale=Float64(row.prop_mass_scale),
        panel_mass_scale=Float64(row.panel_mass_scale),
        bus_dims_scale=Float64(row.bus_dims_scale),
        panel_dims_scale=Float64(row.panel_dims_scale),
        panel_offset_scale=Float64(row.panel_offset_scale),
        srp_cr_scale=Float64(row.srp_cr_scale),
        srp_area_scale=Float64(row.srp_area_scale),
        gravity_variant=String(row.gravity_variant),
        dv_global_scale=Float64(row.dv_global_scale),
        dv_early_scale=Float64(row.dv_early_scale),
        dv_orbit7_bias_mps=Float64(row.dv_orbit7_bias_mps),
        solver_mode=String(row.solver_mode_requested),
        dt_max_orbit_s=Float64(row.dt_max_orbit_requested_s),
        dt_max_atm_s=Float64(row.dt_max_atm_requested_s)
    )
end

function _write_report(
    path::String,
    cfg::TunerConfig,
    quick_df::DataFrame,
    full_df::DataFrame,
    best_row,
    best_manifest::String;
    plot_paths::Vector{String}=String[]
)
    generated = string(now(UTC))
    open(path, "w") do io
        println(io, "# Odyssey Telemetry Tuning Report")
        println(io)
        println(io, "- Generated (UTC): $generated")
        println(io, "- Base manifest: `$(cfg.base_manifest)`")
        println(io, "- Scenario: `$(cfg.scenario_name)`")
        println(io, "- Quick candidates: `$(cfg.quick_candidates)`")
        println(io, "- Full top-k: `$(cfg.full_topk)`")
        println(io, "- Parallel candidates: `$(cfg.parallel_candidates)`")
        println(io, "- Calibration enabled: `$(cfg.calibration_enabled)`")
        println(io, "- Best manifest: `$(best_manifest)`")
        println(io)

        println(io, "## Best Candidate")
        println(io)
        println(io, "- Stage selected: `$(best_row.stage)`")
        println(io, "- Score: `$(best_row.score)`")
        println(io, "- Passes paper limits: `$(best_row.pass)`")
        println(io, "- Peri RMSE [km]: `$(best_row.peri_rmse_km)`")
        println(io, "- Apo RMSE [km]: `$(best_row.apo_rmse_km)`")
        println(io, "- Peri max abs [km]: `$(best_row.peri_max_abs_km)`")
        println(io, "- Apo max abs [km]: `$(best_row.apo_max_abs_km)`")
        println(io)

        println(io, "## Reproduce Best")
        println(io)
        println(io, "```bash")
        println(io, "JULIA_NUM_THREADS=$(Threads.nthreads()) julia --project=.AGORA --startup-file=no scripts/verify_telemetry.jl --profile=full --manifest=$(best_manifest) --enforce=$(cfg.enforce ? "1" : "0")")
        println(io, "```")
        println(io)
        println(io, "## Artifacts")
        println(io)
        println(io, "- Quick results: `$(joinpath(cfg.outdir, "odyssey_tuning_quick_results.csv"))`")
        println(io, "- Full results: `$(joinpath(cfg.outdir, "odyssey_tuning_full_results.csv"))`")
        println(io, "- Combined results: `$(joinpath(cfg.outdir, "odyssey_tuning_all_results.csv"))`")
        println(io, "- Best manifest: `$(best_manifest)`")
        if !isempty(plot_paths)
            println(io)
            println(io, "## Plots")
            println(io)
            for pth in plot_paths
                println(io, "- `$(pth)`")
            end
        end
    end
end

function main_tuner()
    cfg = parse_tuner_cli(ARGS)
    mkpath(cfg.outdir)
    base_manifest = TOML.parsefile(cfg.base_manifest)
    base_scenario = _find_base_scenario(base_manifest, cfg.scenario_name)
    candidates = sample_candidates(cfg)

    println("Odyssey tuner")
    println("outdir=$(cfg.outdir)")
    println("scenario=$(cfg.scenario_name)")
    println("quick_candidates=$(length(candidates)) full_topk=$(cfg.full_topk) parallel_candidates=$(cfg.parallel_candidates) quick_only=$(cfg.quick_only)")

    quick_df = _evaluate_candidates_batch(cfg, :quick, candidates, base_scenario)
    quick_csv = joinpath(cfg.outdir, "odyssey_tuning_quick_results.csv")
    CSV.write(quick_csv, quick_df)

    valid_quick = quick_df[quick_df.success .== true, :]
    if nrow(valid_quick) == 0
        error("All quick candidates failed. Check logs under $(joinpath(cfg.outdir, "logs")).")
    end
    sort!(valid_quick, [:pass, :score, :peri_rmse_km, :apo_rmse_km], rev=[true, false, false, false])

    full_candidates = TuneCandidate[]
    if !cfg.quick_only
        k = min(cfg.full_topk, nrow(valid_quick))
        for i in 1:k
            row = valid_quick[i, :]
            push!(full_candidates, _candidate_from_row(row))
        end
    end
    full_df = _evaluate_candidates_batch(cfg, :full, full_candidates, base_scenario)
    full_csv = joinpath(cfg.outdir, "odyssey_tuning_full_results.csv")
    CSV.write(full_csv, full_df)

    combined_df = vcat(quick_df, full_df; cols=:union)
    combined_csv = joinpath(cfg.outdir, "odyssey_tuning_all_results.csv")
    CSV.write(combined_csv, combined_df)

    best_row = nothing
    if nrow(full_df) > 0
        valid_full = full_df[full_df.success .== true, :]
        if nrow(valid_full) > 0
            sort!(valid_full, [:pass, :score, :peri_rmse_km, :apo_rmse_km], rev=[true, false, false, false])
            best_row = valid_full[1, :]
        else
            sort!(full_df, :score)
            best_row = full_df[1, :]
        end
    else
        best_row = valid_quick[1, :]
    end

    best_cand = _candidate_from_row(best_row)
    best_scenario = _apply_candidate(base_scenario, best_cand; calibration_enabled=cfg.calibration_enabled)
    best_manifest_doc = Dict{String, Any}("version" => 1, "scenarios" => Any[best_scenario])
    best_manifest_path = joinpath(cfg.outdir, "odyssey_tuning_best_manifest.toml")
    open(best_manifest_path, "w") do io
        TOML.print(io, best_manifest_doc)
    end

    plot_paths = _generate_tuning_plots(cfg.outdir, quick_df, full_df)
    report_path = joinpath(cfg.outdir, "odyssey_tuning_report.md")
    _write_report(report_path, cfg, quick_df, full_df, best_row, best_manifest_path; plot_paths=plot_paths)

    println()
    println("Tuning complete.")
    println("Quick results: $quick_csv")
    println("Full results : $full_csv")
    println("Combined CSV : $combined_csv")
    println("Best manifest: $best_manifest_path")
    println("Plots        : $(length(plot_paths))")
    println("Report       : $report_path")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_tuner()
end
