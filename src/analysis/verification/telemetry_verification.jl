module TelemetryVerification

export VerificationRequest, VerificationResult, run_verification, run_verification_cli, run_study

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_DIR = joinpath(REPO_ROOT, "output")
const DEFAULT_MANIFEST_PATH = joinpath(REPO_ROOT, "test", "telemetry_benchmark_manifest.toml")
const STRICT_REL_ORBIT = 1e-7
const STRICT_ABS_ORBIT = 1e-9
const STRICT_REL_ATM = 1e-7
const STRICT_ABS_ATM = 1e-9
const STRICT_DT_ORBIT = 60.0
const STRICT_DT_ATM = 0.2
const TELEMETRY_SOLVER_MAXITERS_QUICK_DEFAULT = 5_000_000
const TELEMETRY_SOLVER_MAXITERS_FULL_DEFAULT = 20_000_000

using Statistics
using Dates
using Printf
using LinearAlgebra
using StaticArrays
using SPICE
using CSV
using DataFrames
using Arrow
using TOML

if isdefined(parentmodule(@__MODULE__), :SimulationEngine)
    using ..SimulationEngine
    const SimulationModel = SimulationEngine.SimulationModel
else
    error(
        "TelemetryVerification requires SimulationEngine in the parent module. " *
        "Load src/simulation/engine/simulation_engine.jl before including telemetry_verification.jl."
    )
end
using .SimulationModel
include(joinpath(REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))
include(joinpath(REPO_ROOT, "src", "examples", "typed_example_utils.jl"))

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EventTolerance = NamedTuple{
    (:max_abs_km, :max_nmae, :max_rmse_km),
    Tuple{Float64, Float64, Float64}
}

Base.@kwdef struct AtmosphereTruthConfig
    assumption_id::String = "gram_default"
    atmosphere_model::String = "GRAM"
    atmosphere_dataset::String = "default"
    space_weather_model::String = "default"
    solar_flux_model::String = "default"
    gram_seed::Int = 1001
    gram_perturbation_scales::NTuple{4, Float64} = (0.0, 0.0, 0.0, 0.0)
    gram_min_relative_step_size::Union{Nothing, Float64} = nothing
    gram_offline_surrogate::String = "off"
    gram_static_grid::Bool = false
    gram_track_cache::Bool = false
    gram_global_lock::String = "on"
    mars_map_year::Union{Nothing, Int} = nothing
    mars_mgcm_dust_levels::Union{Nothing, NTuple{3, Float64}} = nothing
    mars_dust_storm::Union{Nothing, NTuple{6, Float64}} = nothing
    mars_f107::Union{Nothing, Float64} = nothing
    mars_wind_scales::Union{Nothing, NTuple{2, Float64}} = nothing
    mars_mola_heights::Union{Nothing, Bool} = nothing
    mars_min_max::Union{Nothing, Int} = nothing
end

Base.@kwdef struct CalibrationConfig
    enabled::Bool = false
    profiles::Vector{Symbol} = Symbol[:full]
    search_on_quick_subset::Bool = true
    fit_cd_scale::Bool = true
    cd_scale_min::Float64 = 0.85
    cd_scale_max::Float64 = 1.15
    cd_scale_steps::Int = 3
    fit_cr::Bool = true
    cr_min::Float64 = 1.15
    cr_max::Float64 = 1.45
    cr_steps::Int = 3
    fit_bias::Bool = true
    bias_abs_max_km::Float64 = 500.0
    objective::String = "mean_nmae"
end

Base.@kwdef struct SpacecraftConfig
    bus_dims::NTuple{3, Float64}
    panel_dims::NTuple{3, Float64}
    bus_mass_kg::Float64
    panel_mass_each_kg::Float64
    panel_offset_y_m::Float64
    prop_mass_kg::Float64
    id::Int64
end

abstract type AbstractScenarioConfig end

Base.@kwdef struct OrbitEventsScenarioConfig <: AbstractScenarioConfig
    name::String
    planet_name::String
    telemetry_peri_path::String
    telemetry_apo_path::String
    target_orbits_quick::Int
    target_orbits_full::Int
    compare_points_quick::Int
    compare_points_full::Int
    min_eval_points::Int
    units_x::String
    units_y::Dict{String, String}
    tolerances_quick::Dict{String, EventTolerance}
    tolerances_full::Dict{String, EventTolerance}
    initial_time::InitialTime
    ra_m::Float64
    rp_altitude_m::Float64
    i_deg::Float64
    aop_deg::Float64
    raan_deg::Float64
    ta_deg::Float64
    spacecraft::SpacecraftConfig
    gravity_model::Symbol
    gravity_harmonics_degree::Int = 0
    gravity_harmonics_order::Int = 0
    gravity_harmonics_file::String = ""
    nbody_bodies::Vector{String} = String[]
    srp_enabled::Bool = false
    srp_cr::Float64 = 1.3
    srp_area_m2::Float64 = 0.0
    include_wind::Bool = false
    orbit_altitude_mode::Symbol = :vacuum
    maneuver_orbit_numbers::Vector{Int64} = Int64[]
    maneuver_delta_v_mps::Vector{Float64} = Float64[]
    maneuver_thrust_n::Float64 = 0.0
    maneuver_isp_s::Float64 = 0.0
    maneuver_guidance_rate_s::Float64 = 30.0
    maneuver_control_rate_s::Float64 = 10.0
    atmosphere_truth::AtmosphereTruthConfig = AtmosphereTruthConfig()
    calibration::CalibrationConfig = CalibrationConfig()
    EI_km::Float64
end

Base.@kwdef struct TimeAlignedScenarioConfig <: AbstractScenarioConfig
    name::String
    planet_name::String
    telemetry_path::String
    telemetry_time_col::String
    telemetry_altitude_col::String
    telemetry_x_col::String
    telemetry_y_col::String
    telemetry_z_col::String
    telemetry_sma_col::String
    telemetry_ecc_col::String
    telemetry_inc_col::String
    telemetry_aop_col::String
    telemetry_raan_col::String
    telemetry_ta_col::String
    max_points_quick::Int
    max_points_full::Int
    min_eval_points::Int
    units_x::String
    units_y::Dict{String, String}
    tolerances_quick::Dict{String, EventTolerance}
    tolerances_full::Dict{String, EventTolerance}
    initial_time::InitialTime
    spacecraft::SpacecraftConfig
    gravity_model::Symbol
    gravity_harmonics_degree::Int = 0
    gravity_harmonics_order::Int = 0
    gravity_harmonics_file::String = ""
    nbody_bodies::Vector{String} = String[]
    srp_enabled::Bool = false
    srp_cr::Float64 = 1.3
    srp_area_m2::Float64 = 0.0
    include_wind::Bool = false
    orbit_altitude_mode::Symbol = :vacuum
    comparison_mode::Symbol = :time_aligned_state
    extrema_min_separation_s::Float64 = 500.0
    atmosphere_truth::AtmosphereTruthConfig = AtmosphereTruthConfig()
    calibration::CalibrationConfig = CalibrationConfig()
    EI_km::Float64
end

Base.@kwdef struct StudyConfig
    profile::Symbol
    out_summary::String
    out_errors::String
    manifest_path::String
    enforce::Bool
    generate_plots::Bool
end

Base.@kwdef struct VerificationRequest
    profile::Symbol = Symbol(lowercase(get(ENV, "SPACEAGORA_TELEMETRY_PROFILE", "quick")))
    out_summary::String = abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_SUMMARY", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_summary.csv")))
    out_errors::String = abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_ERRORS", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_errors.csv")))
    manifest_path::String = abspath(get(ENV, "SPACEAGORA_TELEMETRY_MANIFEST", DEFAULT_MANIFEST_PATH))
    enforce::Bool = false
    generate_plots::Bool = _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_PLOTS", "1"), true)
end

Base.@kwdef struct VerificationResult
    summary::DataFrame
    errors::DataFrame
    summary_path::String
    errors_path::String
    plots_dir::String
    profile::Symbol
    enforce::Bool
    total_runtime_s::Float64
end

@inline function _safe_parse_bool(raw::AbstractString, default::Bool)::Bool
    token = lowercase(strip(raw))
    if token in ("1", "true", "yes", "on")
        return true
    elseif token in ("0", "false", "no", "off")
        return false
    end
    return default
end

@inline function _parse_positive_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return default
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be a positive integer, got '$raw'"))
    end
    parsed > 0 || throw(ArgumentError("$name must be a positive integer, got $parsed"))
    return parsed
end

@inline function _parse_positive_float_env(name::String)::Union{Nothing, Float64}
    raw = strip(get(ENV, name, ""))
    isempty(raw) && return nothing
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a positive float, got '$raw'"))
    end
    parsed > 0.0 || throw(ArgumentError("$name must be > 0, got $parsed"))
    return parsed
end

@inline function _telemetry_solver_maxiters(profile::Symbol)::Int
    default_val = profile == :quick ? TELEMETRY_SOLVER_MAXITERS_QUICK_DEFAULT : TELEMETRY_SOLVER_MAXITERS_FULL_DEFAULT
    return _parse_positive_int_env("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS", default_val)
end

@inline function _telemetry_solver_retry_maxiters(base_maxiters::Int)::Int
    default_retry = max(base_maxiters * 4, base_maxiters + 1_000_000)
    retry = _parse_positive_int_env("SPACEAGORA_TELEMETRY_SOLVER_MAXITERS_RETRY", default_retry)
    return max(retry, base_maxiters + 1)
end

@inline function _telemetry_solver_mode()::String
    mode = strip(get(ENV, "SPACEAGORA_TELEMETRY_SOLVER_MODE", ""))
    if isempty(mode)
        mode = strip(get(ENV, "SPACEAGORA_SOLVER_MODE", ""))
    end
    return isempty(mode) ? "auto_stiff" : mode
end

@inline _is_maxiters_error(err)::Bool = occursin("MaxIters", sprint(showerror, err))

@inline function _require_key(tbl, key::String, context::String)
    haskey(tbl, key) || throw(ArgumentError("Missing '$key' in $context"))
    return tbl[key]
end

@inline function _require_table(tbl, key::String, context::String)
    value = _require_key(tbl, key, context)
    value isa AbstractDict || throw(ArgumentError("Expected table '$key' in $context"))
    return value
end

@inline function _require_str(tbl, key::String, context::String)::String
    return String(_require_key(tbl, key, context))
end

@inline function _require_float(tbl, key::String, context::String)::Float64
    return Float64(_require_key(tbl, key, context))
end

@inline function _require_int(tbl, key::String, context::String)::Int
    return Int(_require_key(tbl, key, context))
end

@inline function _optional_int(tbl, key::String, default::Int)::Int
    return haskey(tbl, key) ? Int(tbl[key]) : default
end

@inline function _optional_float(tbl, key::String, default::Float64)::Float64
    return haskey(tbl, key) ? Float64(tbl[key]) : default
end

@inline function _optional_bool(tbl, key::String, default::Bool)::Bool
    return haskey(tbl, key) ? Bool(tbl[key]) : default
end

@inline function _optional_str(tbl, key::String, default::String)::String
    return haskey(tbl, key) ? String(tbl[key]) : default
end

@inline function _optional_symbol_vector(tbl, key::String, default::Vector{Symbol})::Vector{Symbol}
    if !haskey(tbl, key)
        return copy(default)
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    out = Symbol[]
    for token in raw
        push!(out, Symbol(lowercase(strip(String(token)))))
    end
    isempty(out) && throw(ArgumentError("Calibration profiles cannot be empty."))
    return out
end

@inline function _require_str_vector(tbl, key::String, context::String)::Vector{String}
    raw = _require_key(tbl, key, context)
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in $context"))
    return [String(v) for v in raw]
end

@inline function _optional_str_vector(tbl, key::String)::Vector{String}
    if !haskey(tbl, key)
        return String[]
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    return [String(v) for v in raw]
end

@inline function _optional_int64_vector(tbl, key::String)::Vector{Int64}
    if !haskey(tbl, key)
        return Int64[]
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    return Int64[Int64(v) for v in raw]
end

@inline function _optional_float64_vector(tbl, key::String)::Vector{Float64}
    if !haskey(tbl, key)
        return Float64[]
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in manifest scenario table."))
    return Float64[Float64(v) for v in raw]
end

@inline function _parse_orbit_altitude_mode(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("vacuum", "spherical", "rp_e")
        return :vacuum
    elseif key in ("true", "oblate", "geodetic")
        return :oblate
    end
    throw(ArgumentError(
        "Unsupported orbit_altitude_mode='$raw' in $context; use vacuum|true."
    ))
end

@inline function _parse_time_aligned_comparison_mode(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("time_aligned_state", "time_aligned", "state")
        return :time_aligned_state
    elseif key in ("orbit_events", "apo_peri", "extrema")
        return :orbit_events
    end
    throw(ArgumentError(
        "Unsupported comparison_mode='$raw' in $context; use time_aligned_state|orbit_events."
    ))
end

function _parse_maneuver_config(tbl, context::String)
    if !haskey(tbl, "maneuvers")
        return (
            orbit_numbers=Int64[],
            delta_v_mps=Float64[],
            thrust_n=0.0,
            isp_s=0.0,
            guidance_rate_s=30.0,
            control_rate_s=10.0
        )
    end
    mtbl = _require_table(tbl, "maneuvers", context)
    orbit_numbers = _optional_int64_vector(mtbl, "orbit_numbers")
    delta_v_mps = _optional_float64_vector(mtbl, "delta_v_mps")
    isempty(orbit_numbers) && throw(ArgumentError("maneuvers.orbit_numbers must be non-empty in $context"))
    length(delta_v_mps) == length(orbit_numbers) || throw(ArgumentError(
        "maneuvers.delta_v_mps length ($(length(delta_v_mps))) must match maneuvers.orbit_numbers length ($(length(orbit_numbers))) in $context"
    ))
    any(v -> v <= 0, orbit_numbers) && throw(ArgumentError("maneuvers.orbit_numbers must be positive integers in $context"))
    return (
        orbit_numbers=orbit_numbers,
        delta_v_mps=delta_v_mps,
        thrust_n=_optional_float(mtbl, "thrust_n", 4.0),
        isp_s=_optional_float(mtbl, "isp_s", 220.0),
        guidance_rate_s=_optional_float(mtbl, "guidance_rate_s", 30.0),
        control_rate_s=_optional_float(mtbl, "control_rate_s", 10.0)
    )
end

@inline function _optional_float_tuple(tbl, key::String, n::Int, context::String)
    if !haskey(tbl, key)
        return nothing
    end
    raw = tbl[key]
    raw isa AbstractVector || throw(ArgumentError("Expected array '$key' in $context"))
    length(raw) == n || throw(ArgumentError("Expected '$key' length=$n in $context, got $(length(raw))"))
    values = ntuple(i -> Float64(raw[i]), n)
    return values
end

function _parse_atmosphere_truth_config(tbl, context::String)::AtmosphereTruthConfig
    if !haskey(tbl, "atmosphere_truth")
        return AtmosphereTruthConfig()
    end
    t = _require_table(tbl, "atmosphere_truth", context)
    assumption_id = _optional_str(t, "assumption_id", "gram_default")
    atmosphere_model = _require_str(t, "atmosphere_model", "$context.atmosphere_truth")
    atmosphere_dataset = _require_str(t, "atmosphere_dataset", "$context.atmosphere_truth")
    space_weather_model = _require_str(t, "space_weather_model", "$context.atmosphere_truth")
    solar_flux_model = _require_str(t, "solar_flux_model", "$context.atmosphere_truth")
    seed = _optional_int(t, "gram_seed", 1001)
    scales_raw = _optional_float_tuple(t, "gram_perturbation_scales", 4, "$context.atmosphere_truth")
    scales = scales_raw === nothing ? (0.0, 0.0, 0.0, 0.0) : (
        scales_raw[1], scales_raw[2], scales_raw[3], scales_raw[4]
    )
    min_step = haskey(t, "gram_min_relative_step_size") ? _optional_float(t, "gram_min_relative_step_size", 0.0) : nothing
    offline_surrogate = lowercase(_optional_str(t, "gram_offline_surrogate", "off"))
    offline_surrogate in ("off", "on", "auto") || throw(ArgumentError(
        "Unsupported atmosphere_truth.gram_offline_surrogate='$offline_surrogate' in $context; use off|on|auto."
    ))
    global_lock = lowercase(_optional_str(t, "gram_global_lock", "on"))
    global_lock in ("on", "off") || throw(ArgumentError(
        "Unsupported atmosphere_truth.gram_global_lock='$global_lock' in $context; use on|off."
    ))
    mars_mgcm_raw = _optional_float_tuple(t, "mars_mgcm_dust_levels", 3, "$context.atmosphere_truth")
    mars_dust_raw = _optional_float_tuple(t, "mars_dust_storm", 6, "$context.atmosphere_truth")
    mars_wind_raw = _optional_float_tuple(t, "mars_wind_scales", 2, "$context.atmosphere_truth")
    return AtmosphereTruthConfig(
        assumption_id=assumption_id,
        atmosphere_model=atmosphere_model,
        atmosphere_dataset=atmosphere_dataset,
        space_weather_model=space_weather_model,
        solar_flux_model=solar_flux_model,
        gram_seed=seed,
        gram_perturbation_scales=scales,
        gram_min_relative_step_size=min_step,
        gram_offline_surrogate=offline_surrogate,
        gram_static_grid=_optional_bool(t, "gram_static_grid", false),
        gram_track_cache=_optional_bool(t, "gram_track_cache", false),
        gram_global_lock=global_lock,
        mars_map_year=haskey(t, "mars_map_year") ? _optional_int(t, "mars_map_year", 0) : nothing,
        mars_mgcm_dust_levels=mars_mgcm_raw === nothing ? nothing : (
            mars_mgcm_raw[1], mars_mgcm_raw[2], mars_mgcm_raw[3]
        ),
        mars_dust_storm=mars_dust_raw === nothing ? nothing : (
            mars_dust_raw[1], mars_dust_raw[2], mars_dust_raw[3], mars_dust_raw[4], mars_dust_raw[5], mars_dust_raw[6]
        ),
        mars_f107=haskey(t, "mars_f107") ? _optional_float(t, "mars_f107", 0.0) : nothing,
        mars_wind_scales=mars_wind_raw === nothing ? nothing : (mars_wind_raw[1], mars_wind_raw[2]),
        mars_mola_heights=haskey(t, "mars_mola_heights") ? _optional_bool(t, "mars_mola_heights", true) : nothing,
        mars_min_max=haskey(t, "mars_min_max") ? _optional_int(t, "mars_min_max", 0) : nothing
    )
end

function _parse_calibration_config(tbl, context::String)::CalibrationConfig
    if !haskey(tbl, "calibration")
        return CalibrationConfig()
    end
    c = _require_table(tbl, "calibration", context)
    profiles = _optional_symbol_vector(c, "profiles", Symbol[:full])
    objective = lowercase(_optional_str(c, "objective", "mean_nmae"))
    objective in ("mean_nmae", "mean_rmse_km", "max_nmae") || throw(ArgumentError(
        "Unsupported calibration.objective='$objective' in $context; use mean_nmae|mean_rmse_km|max_nmae."
    ))
    cd_steps = _optional_int(c, "cd_scale_steps", 3)
    cr_steps = _optional_int(c, "cr_steps", 3)
    cd_steps > 0 || throw(ArgumentError("calibration.cd_scale_steps must be > 0 in $context"))
    cr_steps > 0 || throw(ArgumentError("calibration.cr_steps must be > 0 in $context"))
    cd_min = _optional_float(c, "cd_scale_min", 0.85)
    cd_max = _optional_float(c, "cd_scale_max", 1.15)
    cr_min = _optional_float(c, "cr_min", 1.15)
    cr_max = _optional_float(c, "cr_max", 1.45)
    cd_max >= cd_min || throw(ArgumentError("calibration.cd_scale_max must be >= cd_scale_min in $context"))
    cr_max >= cr_min || throw(ArgumentError("calibration.cr_max must be >= cr_min in $context"))
    return CalibrationConfig(
        enabled=_optional_bool(c, "enabled", false),
        profiles=profiles,
        search_on_quick_subset=_optional_bool(c, "search_on_quick_subset", true),
        fit_cd_scale=_optional_bool(c, "fit_cd_scale", true),
        cd_scale_min=cd_min,
        cd_scale_max=cd_max,
        cd_scale_steps=cd_steps,
        fit_cr=_optional_bool(c, "fit_cr", true),
        cr_min=cr_min,
        cr_max=cr_max,
        cr_steps=cr_steps,
        fit_bias=_optional_bool(c, "fit_bias", true),
        bias_abs_max_km=_optional_float(c, "bias_abs_max_km", 500.0),
        objective=objective
    )
end

@inline function _resolve_repo_path(path::String)::String
    return isabspath(path) ? path : normpath(joinpath(REPO_ROOT, path))
end

@inline function _parse_gravity_model(raw::String, context::String)::Symbol
    key = lowercase(strip(raw))
    if key in ("inverse_squared", "is2", "point_mass")
        return :inverse_squared
    elseif key in ("inverse_squared_j2", "is2_j2", "j2")
        return :inverse_squared_j2
    end
    throw(ArgumentError("Unsupported gravity_model='$raw' in $context"))
end

function _parse_initial_time(tbl, context::String)::InitialTime
    return InitialTime(
        year=_require_int(tbl, "year", context),
        month=_require_int(tbl, "month", context),
        day=_require_int(tbl, "day", context),
        hour=_require_int(tbl, "hour", context),
        minute=_require_int(tbl, "minute", context),
        second=_require_float(tbl, "second", context)
    )
end

function _parse_vec3(tbl, key::String, context::String)::NTuple{3, Float64}
    raw = _require_key(tbl, key, context)
    raw isa AbstractVector || throw(ArgumentError("Expected 3-element array '$key' in $context"))
    length(raw) == 3 || throw(ArgumentError("Expected 3 values for '$key' in $context, got $(length(raw))"))
    return (Float64(raw[1]), Float64(raw[2]), Float64(raw[3]))
end

function _parse_spacecraft_config(tbl, context::String)::SpacecraftConfig
    stbl = _require_table(tbl, "spacecraft", context)
    return SpacecraftConfig(
        bus_dims=_parse_vec3(stbl, "bus_dims_m", "$context.spacecraft"),
        panel_dims=_parse_vec3(stbl, "panel_dims_m", "$context.spacecraft"),
        bus_mass_kg=_require_float(stbl, "bus_mass_kg", "$context.spacecraft"),
        panel_mass_each_kg=_require_float(stbl, "panel_mass_each_kg", "$context.spacecraft"),
        panel_offset_y_m=_require_float(stbl, "panel_offset_y_m", "$context.spacecraft"),
        prop_mass_kg=_require_float(stbl, "prop_mass_kg", "$context.spacecraft"),
        id=Int64(_require_int(stbl, "id", "$context.spacecraft"))
    )
end

function _parse_units(tbl, events::Vector{String}, context::String)::Tuple{String, Dict{String, String}}
    utbl = _require_table(tbl, "units", context)
    x_units = _require_str(utbl, "x", "$context.units")
    y_units = Dict{String, String}()
    for event in events
        y_units[event] = _require_str(utbl, event, "$context.units")
    end
    return x_units, y_units
end

function _parse_tolerances(tbl, key::String, events::Vector{String}, context::String)::Dict{String, EventTolerance}
    ttbl = _require_table(tbl, key, context)
    out = Dict{String, EventTolerance}()
    for event in events
        etbl = _require_table(ttbl, event, "$context.$key")
        out[event] = (
            max_abs_km=_require_float(etbl, "max_abs_km", "$context.$key.$event"),
            max_nmae=_require_float(etbl, "max_nmae", "$context.$key.$event"),
            max_rmse_km=_optional_float(etbl, "max_rmse_km", Inf)
        )
    end
    return out
end

function _load_scenarios_from_manifest(manifest_path::String)::Vector{AbstractScenarioConfig}
    doc = TOML.parsefile(manifest_path)
    raw = _require_key(doc, "scenarios", "manifest")
    raw isa AbstractVector || throw(ArgumentError("Manifest field 'scenarios' must be an array of tables."))

    scenarios = AbstractScenarioConfig[]
    for (idx, entry) in pairs(raw)
        entry isa AbstractDict || throw(ArgumentError("manifest.scenarios[$idx] must be a table."))
        tbl = entry
        context = "manifest.scenarios[$idx]"

        name = _require_str(tbl, "name", context)
        kind = lowercase(_require_str(tbl, "kind", context))
        planet_name = lowercase(_require_str(tbl, "planet", context))
        events = _require_str_vector(tbl, "events", context)
        units_x, units_y = _parse_units(tbl, events, context)
        tolerances_quick = _parse_tolerances(tbl, "tolerances_quick", events, context)
        tolerances_full = _parse_tolerances(tbl, "tolerances_full", events, context)
        min_eval_points = _require_int(tbl, "min_eval_points", context)
        initial_time = _parse_initial_time(_require_table(tbl, "initial_time", context), "$context.initial_time")
        spacecraft = _parse_spacecraft_config(tbl, context)
        gravity_model = _parse_gravity_model(_require_str(tbl, "gravity_model", context), context)
        EI_km = _require_float(tbl, "EI_km", context)
        gravity_harmonics_degree = _optional_int(tbl, "gravity_harmonics_degree", 0)
        gravity_harmonics_order = _optional_int(tbl, "gravity_harmonics_order", 0)
        raw_harmonics_file = _optional_str(tbl, "gravity_harmonics_file", "")
        gravity_harmonics_file = isempty(strip(raw_harmonics_file)) ? "" : _resolve_repo_path(raw_harmonics_file)
        nbody_bodies = _optional_str_vector(tbl, "nbody_bodies")
        srp_enabled = _optional_bool(tbl, "srp_enabled", false)
        srp_cr = _optional_float(tbl, "srp_cr", 1.3)
        srp_area_m2 = _optional_float(tbl, "srp_area_m2", 0.0)
        include_wind = _optional_bool(tbl, "include_wind", false)
        orbit_altitude_mode = _parse_orbit_altitude_mode(_optional_str(tbl, "orbit_altitude_mode", "vacuum"), context)
        maneuver = _parse_maneuver_config(tbl, context)
        atmosphere_truth = _parse_atmosphere_truth_config(tbl, context)
        calibration = _parse_calibration_config(tbl, context)

        if kind == "orbit_events"
            push!(scenarios, OrbitEventsScenarioConfig(
                name=name,
                planet_name=planet_name,
                telemetry_peri_path=_resolve_repo_path(_require_str(tbl, "telemetry_peri", context)),
                telemetry_apo_path=_resolve_repo_path(_require_str(tbl, "telemetry_apo", context)),
                target_orbits_quick=_require_int(tbl, "target_orbits_quick", context),
                target_orbits_full=_require_int(tbl, "target_orbits_full", context),
                compare_points_quick=_require_int(tbl, "compare_points_quick", context),
                compare_points_full=_require_int(tbl, "compare_points_full", context),
                min_eval_points=min_eval_points,
                units_x=units_x,
                units_y=units_y,
                tolerances_quick=tolerances_quick,
                tolerances_full=tolerances_full,
                initial_time=initial_time,
                ra_m=_require_float(tbl, "ra_m", context),
                rp_altitude_m=_require_float(tbl, "rp_altitude_m", context),
                i_deg=_require_float(tbl, "i_deg", context),
                aop_deg=_require_float(tbl, "aop_deg", context),
                raan_deg=_require_float(tbl, "raan_deg", context),
                ta_deg=_require_float(tbl, "ta_deg", context),
                spacecraft=spacecraft,
                gravity_model=gravity_model,
                gravity_harmonics_degree=gravity_harmonics_degree,
                gravity_harmonics_order=gravity_harmonics_order,
                gravity_harmonics_file=gravity_harmonics_file,
                nbody_bodies=nbody_bodies,
                srp_enabled=srp_enabled,
                srp_cr=srp_cr,
                srp_area_m2=srp_area_m2,
                include_wind=include_wind,
                orbit_altitude_mode=orbit_altitude_mode,
                maneuver_orbit_numbers=maneuver.orbit_numbers,
                maneuver_delta_v_mps=maneuver.delta_v_mps,
                maneuver_thrust_n=maneuver.thrust_n,
                maneuver_isp_s=maneuver.isp_s,
                maneuver_guidance_rate_s=maneuver.guidance_rate_s,
                maneuver_control_rate_s=maneuver.control_rate_s,
                atmosphere_truth=atmosphere_truth,
                calibration=calibration,
                EI_km=EI_km
            ))
        elseif kind == "time_aligned_state"
            ctbl = _require_table(tbl, "telemetry_columns", context)
            comparison_mode = _parse_time_aligned_comparison_mode(
                _optional_str(tbl, "comparison_mode", "time_aligned_state"),
                context
            )
            extrema_min_separation_s = _optional_float(tbl, "extrema_min_separation_s", 500.0)
            extrema_min_separation_s > 0.0 || throw(ArgumentError(
                "extrema_min_separation_s must be > 0.0 in $context"
            ))
            push!(scenarios, TimeAlignedScenarioConfig(
                name=name,
                planet_name=planet_name,
                telemetry_path=_resolve_repo_path(_require_str(tbl, "telemetry", context)),
                telemetry_time_col=_require_str(ctbl, "time", "$context.telemetry_columns"),
                telemetry_altitude_col=_require_str(ctbl, "altitude", "$context.telemetry_columns"),
                telemetry_x_col=_require_str(ctbl, "x", "$context.telemetry_columns"),
                telemetry_y_col=_require_str(ctbl, "y", "$context.telemetry_columns"),
                telemetry_z_col=_require_str(ctbl, "z", "$context.telemetry_columns"),
                telemetry_sma_col=_require_str(ctbl, "sma", "$context.telemetry_columns"),
                telemetry_ecc_col=_require_str(ctbl, "ecc", "$context.telemetry_columns"),
                telemetry_inc_col=_require_str(ctbl, "inc", "$context.telemetry_columns"),
                telemetry_aop_col=_require_str(ctbl, "aop", "$context.telemetry_columns"),
                telemetry_raan_col=_require_str(ctbl, "raan", "$context.telemetry_columns"),
                telemetry_ta_col=_require_str(ctbl, "ta", "$context.telemetry_columns"),
                max_points_quick=_require_int(tbl, "max_points_quick", context),
                max_points_full=_require_int(tbl, "max_points_full", context),
                min_eval_points=min_eval_points,
                units_x=units_x,
                units_y=units_y,
                tolerances_quick=tolerances_quick,
                tolerances_full=tolerances_full,
                initial_time=initial_time,
                spacecraft=spacecraft,
                gravity_model=gravity_model,
                gravity_harmonics_degree=gravity_harmonics_degree,
                gravity_harmonics_order=gravity_harmonics_order,
                gravity_harmonics_file=gravity_harmonics_file,
                nbody_bodies=nbody_bodies,
                srp_enabled=srp_enabled,
                srp_cr=srp_cr,
                srp_area_m2=srp_area_m2,
                include_wind=include_wind,
                orbit_altitude_mode=orbit_altitude_mode,
                comparison_mode=comparison_mode,
                extrema_min_separation_s=extrema_min_separation_s,
                atmosphere_truth=atmosphere_truth,
                calibration=calibration,
                EI_km=EI_km
            ))
        else
            throw(ArgumentError("Unsupported scenario kind '$kind' in $context"))
        end
    end

    isempty(scenarios) && throw(ArgumentError("No scenarios defined in manifest: $manifest_path"))
    return scenarios
end

function parse_cli(args::Vector{String})::StudyConfig
    profile = Symbol(lowercase(get(ENV, "SPACEAGORA_TELEMETRY_PROFILE", "quick")))
    out_summary = get(ENV, "SPACEAGORA_TELEMETRY_OUT_SUMMARY", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_summary.csv"))
    out_errors = get(ENV, "SPACEAGORA_TELEMETRY_OUT_ERRORS", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_errors.csv"))
    manifest_path = get(ENV, "SPACEAGORA_TELEMETRY_MANIFEST", DEFAULT_MANIFEST_PATH)
    enforce = _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_ENFORCE", "0"), false)
    generate_plots = _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_PLOTS", "1"), true)

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("quick", "full")
            profile = Symbol(arg)
        elseif startswith(arg, "--profile=")
            profile = Symbol(lowercase(split(arg, "=", limit=2)[2]))
        elseif startswith(arg, "--out-summary=")
            out_summary = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--out-errors=")
            out_errors = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--manifest=")
            manifest_path = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--enforce=")
            enforce = _safe_parse_bool(split(arg, "=", limit=2)[2], enforce)
        elseif startswith(arg, "--plots=")
            generate_plots = _safe_parse_bool(split(arg, "=", limit=2)[2], generate_plots)
        else
            throw(ArgumentError("Unsupported argument '$arg'. Supported: [quick|full], --profile=..., --manifest=..., --out-summary=..., --out-errors=..., --enforce=true|false, --plots=true|false"))
        end
        i += 1
    end

    profile in (:quick, :full) || throw(ArgumentError("Unsupported profile '$profile'. Use quick or full."))
    return StudyConfig(
        profile=profile,
        out_summary=abspath(out_summary),
        out_errors=abspath(out_errors),
        manifest_path=abspath(manifest_path),
        enforce=enforce,
        generate_plots=generate_plots
    )
end

@inline function _study_config(request::VerificationRequest)::StudyConfig
    profile = Symbol(lowercase(String(request.profile)))
    profile in (:quick, :full) || throw(ArgumentError("Unsupported request.profile '$profile'. Use quick or full."))
    return StudyConfig(
        profile=profile,
        out_summary=abspath(request.out_summary),
        out_errors=abspath(request.out_errors),
        manifest_path=abspath(request.manifest_path),
        enforce=request.enforce,
        generate_plots=request.generate_plots
    )
end

@inline function _request_from_study_config(cfg::StudyConfig)::VerificationRequest
    return VerificationRequest(
        profile=cfg.profile,
        out_summary=cfg.out_summary,
        out_errors=cfg.out_errors,
        manifest_path=cfg.manifest_path,
        enforce=cfg.enforce,
        generate_plots=cfg.generate_plots
    )
end

@inline function _period_seconds(planet, ra::Float64, rp::Float64)::Float64
    a = 0.5 * (ra + rp)
    return 2π * sqrt(a^3 / planet.μ)
end

@inline function _planet_from_name(planet_name::String)
    key = lowercase(strip(planet_name))
    if key == "mars"
        return Mars("", SPICE_PATH)
    elseif key == "venus"
        return Venus("", SPICE_PATH)
    elseif key == "earth"
        return Earth("", SPICE_PATH)
    end
    throw(ArgumentError("Unsupported planet '$planet_name' in telemetry benchmark."))
end

@inline function _base_gravity_effector(gravity_model::Symbol)
    if gravity_model == :inverse_squared
        return InverseSquaredGravityModel()
    elseif gravity_model == :inverse_squared_j2
        return InverseSquaredJ2GravityModel()
    end
    throw(ArgumentError("Unsupported gravity model '$gravity_model'"))
end

@inline function _harmonics_order(degree::Int, order::Int)::Int
    if degree <= 0
        return 0
    end
    if order <= 0
        return degree
    end
    return min(order, degree)
end

@inline function _nbody_primary_name(planet_name::String)::String
    key = lowercase(strip(planet_name))
    if key == "earth"
        return "Earth"
    elseif key == "mars"
        return "Mars"
    elseif key == "venus"
        return "Venus"
    elseif key == "titan"
        return "Titan"
    end
    throw(ArgumentError("Unsupported N-body primary planet '$planet_name'"))
end

struct ScaledAerodynamicCoefficientfM <: AbstractForceTorqueModel
    model::AerodynamicCoefficientfM
    cd_scale::Float64
    function ScaledAerodynamicCoefficientfM(model::AerodynamicCoefficientfM, cd_scale::Float64)
        cd_scale > 0.0 || throw(ArgumentError("ScaledAerodynamicCoefficientfM.cd_scale must be > 0.0, got $cd_scale"))
        return new(model, cd_scale)
    end
end

@inline _dynamic_effector_threadsafe(::ScaledAerodynamicCoefficientfM)::Bool = true

function SimulationModel.calcForceTorque(
    model::ScaledAerodynamicCoefficientfM,
    x::AbstractVector,
    param::SimulationModel.ODEParams,
    i::Int64
)
    f, τ = SimulationModel.calcForceTorque(model.model, x, param, i)
    return model.cd_scale .* f, model.cd_scale .* τ
end

function _scenario_dynamic_effectors(
    cfg::AbstractScenarioConfig,
    planet,
    spacecraft;
    cd_scale::Float64=1.0,
    cr_override::Union{Nothing, Float64}=nothing
)
    effectors = Any[_base_gravity_effector(cfg.gravity_model)]

    if cfg.gravity_harmonics_degree > 0
        harmonics_file = cfg.gravity_harmonics_file
        isempty(harmonics_file) && throw(ArgumentError(
            "Scenario $(cfg.name) sets gravity_harmonics_degree=$(cfg.gravity_harmonics_degree) but does not provide gravity_harmonics_file."
        ))
        isfile(harmonics_file) || throw(ArgumentError("Missing gravity_harmonics_file for $(cfg.name): $harmonics_file"))
        push!(
            effectors,
            GravitationalHarmonicsModel(
                cfg.gravity_harmonics_degree,
                _harmonics_order(cfg.gravity_harmonics_degree, cfg.gravity_harmonics_order),
                harmonics_file,
                planet
            )
        )
    end

    if !isempty(cfg.nbody_bodies)
        push!(
            effectors,
            NBodyGravityModel(
                body_names=Tuple(cfg.nbody_bodies),
                primary_body_name=_nbody_primary_name(cfg.planet_name),
                planet=planet
            )
        )
    end

    if cfg.srp_enabled
        area_m2 = cfg.srp_area_m2 > 0.0 ? cfg.srp_area_m2 : Float64(spacecraft.root.ref_area)
        area_m2 > 0.0 || throw(ArgumentError("Scenario $(cfg.name) SRP area must be > 0.0 m^2"))
        cr_value = cr_override === nothing ? cfg.srp_cr : Float64(cr_override)
        push!(effectors, SolarRadiationPressureModel(cr_value, area_m2))
    end

    aero = AerodynamicCoefficientfM()
    if isapprox(cd_scale, 1.0; rtol=0.0, atol=1e-12)
        push!(effectors, aero)
    else
        push!(effectors, ScaledAerodynamicCoefficientfM(aero, cd_scale))
    end
    return Tuple(effectors)
end

Base.@kwdef struct _GRAMOfflineSurrogateFallbackBase
    planet_name::String
end

function SimulationModel.EnvironmentModels._gram_point_density(
    model::_GRAMOfflineSurrogateFallbackBase,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    # Library-less fallback cannot query native GRAM point model; return vacuum-like fallback.
    return 0.0, 200.0, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function _is_gram_library_missing_error(err)::Bool
    msg = lowercase(sprint(showerror, err))
    return occursin("gram shared library", msg) || occursin("gram_lib", msg)
end

@inline function _libraryless_gram_surrogate_enabled()::Bool
    return _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_ALLOW_GRAM_OFFLINE_NO_LIB", "1"), true)
end

function _try_libraryless_gram_surrogate(planet_name::String)
    _libraryless_gram_surrogate_enabled() || return nothing
    planet_key = lowercase(strip(planet_name))
    surrogate_file = try
        Base.invokelatest(SimulationModel.EnvironmentModels._gram_default_surrogate_file, planet_key)
    catch
        ""
    end
    isempty(surrogate_file) && return nothing
    isfile(surrogate_file) || return nothing
    base_model = _GRAMOfflineSurrogateFallbackBase(planet_name=planet_key)
    return GRAMAtmosphereModelSurrogate(base_model, surrogate_file, nothing)
end

function _make_required_gram_density_model(
    planet_name::String,
    initial_time::InitialTime,
    truth::AtmosphereTruthConfig
)
    try
        return Base.invokelatest(
            GRAMAtmosphereModel;
            planet_name=planet_name,
            initial_time=initial_time,
            seed=truth.gram_seed,
            gram_min_relative_step_size=truth.gram_min_relative_step_size,
            gram_perturbation_scales=truth.gram_perturbation_scales,
            mars_map_year=truth.mars_map_year,
            mars_mgcm_dust_levels=truth.mars_mgcm_dust_levels,
            mars_dust_storm=truth.mars_dust_storm,
            mars_f107=truth.mars_f107,
            mars_wind_scales=truth.mars_wind_scales,
            mars_mola_heights=truth.mars_mola_heights,
            mars_min_max=truth.mars_min_max
        )
    catch err
        msg = sprint(showerror, err)
        if _is_gram_library_missing_error(err)
            offline_model = _try_libraryless_gram_surrogate(planet_name)
            if offline_model !== nothing
                @warn "GRAM shared library unavailable; using library-less GRAM offline surrogate fallback for telemetry." planet=planet_name surrogate_file=offline_model.surrogate_file
                return offline_model
            end
        end
        throw(ErrorException("GRAM atmosphere initialization failed for planet=$planet_name initial_time=$(initial_time): $msg"))
    end
end

@inline function _make_spacecraft(cfg::SpacecraftConfig, ic::InitialCondition)
    return make_three_body_spacecraft(
        bus_dims=cfg.bus_dims,
        panel_dims=cfg.panel_dims,
        bus_mass=cfg.bus_mass_kg,
        panel_mass_each=cfg.panel_mass_each_kg,
        panel_offset_y=cfg.panel_offset_y_m,
        ic=ic,
        prop_mass=cfg.prop_mass_kg,
        id=cfg.id
    )
end

function _with_environment_wind(args::SimulationConfiguration, include_wind::Bool)::SimulationConfiguration
    env = args.environment_model
    env_updated = EnvironmentModel(
        planet=env.planet,
        EI=env.EI,
        density_model=env.density_model,
        thermal_model=env.thermal_model,
        topography=env.topography,
        topo_degree=env.topo_degree,
        topo_order=env.topo_order,
        wind=include_wind
    )
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration,
        environment_model=env_updated,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

@inline function _has_campaign_maneuvers(cfg::OrbitEventsScenarioConfig)::Bool
    return !isempty(cfg.maneuver_orbit_numbers)
end

function _with_campaign_maneuvers(args::SimulationConfiguration, cfg::OrbitEventsScenarioConfig)::SimulationConfiguration
    !_has_campaign_maneuvers(cfg) && return args

    n_sats = length(args.dynamics_model.spacecraft)
    n_sats == 1 || throw(ArgumentError(
        "Telemetry campaign maneuvers currently support exactly one spacecraft; got $n_sats for scenario $(cfg.name)."
    ))
    cfg.maneuver_thrust_n > 0.0 || throw(ArgumentError("maneuvers.thrust_n must be > 0.0 for scenario $(cfg.name)"))
    cfg.maneuver_isp_s > 0.0 || throw(ArgumentError("maneuvers.isp_s must be > 0.0 for scenario $(cfg.name)"))
    cfg.maneuver_guidance_rate_s > 0.0 || throw(ArgumentError("maneuvers.guidance_rate_s must be > 0.0 for scenario $(cfg.name)"))
    cfg.maneuver_control_rate_s > 0.0 || throw(ArgumentError("maneuvers.control_rate_s must be > 0.0 for scenario $(cfg.name)"))

    thruster = BaseThrusterModel(
        thrust=fill(cfg.maneuver_thrust_n, n_sats),
        direction=fill(0.0, n_sats),
        Δv=fill(0.0, n_sats),
        start_burn_time=fill(-1.0, n_sats),
        stop_burn_time=fill(-1.0, n_sats),
        Isp=fill(cfg.maneuver_isp_s, n_sats)
    )
    guidance_effector = AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=cfg.maneuver_orbit_numbers,
        maneuver_Δv=cfg.maneuver_delta_v_mps
    )
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=args.mission_configuration,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=GuidanceModel(
            guidance_effectors=(guidance_effector,),
            guidance_rates=[cfg.maneuver_guidance_rate_s]
        ),
        navigation_model=args.navigation_model,
        control_model=ControlModel(
            control_effectors=(thruster,),
            control_rates=[cfg.maneuver_control_rate_s]
        ),
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

function _with_orbit_mission(
    args::SimulationConfiguration,
    target_orbits::Int,
    mission_time_s::Float64
)::SimulationConfiguration
    mc = args.mission_configuration
    mission_cfg = MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=mc.keplerian,
        number_of_orbits=target_orbits,
        mission_time=mission_time_s,
        orientation_sim=mc.orientation_sim,
        num_steps_to_save=mc.num_steps_to_save
    )
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=args.simulation_settings,
        mission_configuration=mission_cfg,
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=args.integration_tolerances
    )
end

function _make_orbit_args(
    cfg::OrbitEventsScenarioConfig,
    target_orbits::Int;
    cd_scale::Float64=1.0,
    cr_override::Union{Nothing, Float64}=nothing
)::SimulationConfiguration
    planet = _planet_from_name(cfg.planet_name)
    rp_m = planet.Rp_e + cfg.rp_altitude_m
    ic = InitialCondition(
        ra=cfg.ra_m,
        rp=rp_m,
        i=cfg.i_deg,
        ω=cfg.aop_deg,
        Ω=cfg.raan_deg,
        ν=cfg.ta_deg
    )

    spacecraft = _make_spacecraft(cfg.spacecraft, ic)
    dynamic_effectors = _scenario_dynamic_effectors(
        cfg,
        planet,
        spacecraft;
        cd_scale=cd_scale,
        cr_override=cr_override
    )
    density_model = _make_required_gram_density_model(cfg.planet_name, cfg.initial_time, cfg.atmosphere_truth)

    period_s = _period_seconds(planet, cfg.ra_m, rp_m)
    mission_time_s = target_orbits * period_s

    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time_s,
        initial_time=cfg.initial_time,
        dynamic_effectors=dynamic_effectors,
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=cfg.EI_km,
        verbose=false
    )
    with_orbits = _with_orbit_mission(base, target_orbits, mission_time_s)
    with_wind = _with_environment_wind(with_orbits, cfg.include_wind)
    return _with_campaign_maneuvers(with_wind, cfg)
end

function _make_time_aligned_args(
    cfg::TimeAlignedScenarioConfig,
    mission_time_s::Float64,
    ic::InitialCondition;
    cd_scale::Float64=1.0,
    cr_override::Union{Nothing, Float64}=nothing
)::SimulationConfiguration
    planet = _planet_from_name(cfg.planet_name)
    spacecraft = _make_spacecraft(cfg.spacecraft, ic)
    dynamic_effectors = _scenario_dynamic_effectors(
        cfg,
        planet,
        spacecraft;
        cd_scale=cd_scale,
        cr_override=cr_override
    )
    density_model = _make_required_gram_density_model(cfg.planet_name, cfg.initial_time, cfg.atmosphere_truth)

    base = make_example_config(
        planet=planet,
        spacecraft=spacecraft,
        mission_time=mission_time_s,
        initial_time=cfg.initial_time,
        dynamic_effectors=dynamic_effectors,
        density_model=density_model,
        orientation_sim=false,
        keplerian=false,
        EI_km=cfg.EI_km,
        verbose=false
    )
    return _with_environment_wind(base, cfg.include_wind)
end

@inline function _has_high_fidelity_effectors(args::SimulationConfiguration)::Bool
    return any(
        eff -> (
            eff isa GravitationalHarmonicsModel ||
            eff isa NBodyGravityModel ||
            eff isa SolarRadiationPressureModel
        ),
        args.dynamics_model.dynamic_effectors
    )
end

function _with_study_settings(args::SimulationConfiguration; quick::Bool=false)::SimulationConfiguration
    hf = _has_high_fidelity_effectors(args)
    rel_orbit = min(quick ? 5e-7 : 1e-7, STRICT_REL_ORBIT)
    abs_orbit = min(quick ? 5e-9 : 1e-9, STRICT_ABS_ORBIT)
    rel_atm = min(quick ? 1e-6 : 1e-7, STRICT_REL_ATM)
    abs_atm = min(quick ? 1e-8 : 1e-9, STRICT_ABS_ATM)
    dt_orbit_base = min(quick ? (hf ? 180.0 : 240.0) : (hf ? 60.0 : 120.0), STRICT_DT_ORBIT)
    dt_atm_base = min(quick ? (hf ? 2.0 : 5.0) : (hf ? 0.2 : 0.5), STRICT_DT_ATM)
    dt_orbit_env = _parse_positive_float_env("SPACEAGORA_TELEMETRY_DT_MAX_ORBIT")
    dt_atm_env = _parse_positive_float_env("SPACEAGORA_TELEMETRY_DT_MAX_ATM")
    dt_orbit = dt_orbit_env === nothing ? dt_orbit_base : min(dt_orbit_env, STRICT_DT_ORBIT)
    dt_atm = dt_atm_env === nothing ? dt_atm_base : min(dt_atm_env, STRICT_DT_ATM)
    return SimulationConfiguration(
        file_paths=args.file_paths,
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            results_directory=args.simulation_settings.results_directory,
            generate_plots=false,
            generate_filenames=false,
            normalize=false,
            save_csv=true
        ),
        mission_configuration=MissionConfiguration(
            mission_type=args.mission_configuration.mission_type,
            keplerian=args.mission_configuration.keplerian,
            number_of_orbits=args.mission_configuration.number_of_orbits,
            mission_time=args.mission_configuration.mission_time,
            orientation_sim=args.mission_configuration.orientation_sim,
            num_steps_to_save=2000
        ),
        environment_model=args.environment_model,
        dynamics_model=args.dynamics_model,
        guidance_model=args.guidance_model,
        navigation_model=args.navigation_model,
        control_model=args.control_model,
        initial_time=args.initial_time,
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=rel_orbit,
            abstol_orbit=abs_orbit,
            dt_max_orbit=dt_orbit,
            reltol_atmosphere=rel_atm,
            abstol_atmosphere=abs_atm,
            dt_max_atmosphere=dt_atm
        )
    )
end

function _save_fields_for_study()
    pos_getter = (u, t, integrator) -> begin
        out = Vector{SVector{3, Float64}}(undef, length(u.sc))
        @inbounds for i in eachindex(u.sc)
            out[i] = SVector{3, Float64}(u.sc[i].pos)
        end
        return out
    end
    vel_getter = (u, t, integrator) -> begin
        out = Vector{SVector{3, Float64}}(undef, length(u.sc))
        @inbounds for i in eachindex(u.sc)
            out[i] = SVector{3, Float64}(u.sc[i].vel)
        end
        return out
    end

    return [
        SaveField(:position, pos_getter; per_satellite=true, column_prefix="pos"),
        SaveField(:velocity, vel_getter; per_satellite=true, column_prefix="vel")
    ]
end

@inline function _to_float_vector(values, context::String)::Vector{Float64}
    out = Vector{Float64}(undef, length(values))
    @inbounds for i in eachindex(values)
        v = values[i]
        v === missing && throw(ArgumentError("Missing value in $context at index $i"))
        out[i] = Float64(v)
    end
    return out
end

@inline function _require_column(df::DataFrame, candidates::Vector{String}, context::String)::Vector{Float64}
    for col in candidates
        if col in names(df)
            return _to_float_vector(df[!, col], "$context:$col")
        end
    end
    throw(ArgumentError("Missing required column for $context. Tried: $(join(candidates, ", "))"))
end

function _extract_extrema_series(df::DataFrame, planet, altitude_mode::Symbol)
    x = _require_column(df, ["sc1_pos_1", "sc1_position_1"], "position-x")
    y = _require_column(df, ["sc1_pos_2", "sc1_position_2"], "position-y")
    z = _require_column(df, ["sc1_pos_3", "sc1_position_3"], "position-z")
    vx = _require_column(df, ["sc1_vel_1", "sc1_velocity_1"], "velocity-x")
    vy = _require_column(df, ["sc1_vel_2", "sc1_velocity_2"], "velocity-y")
    vz = _require_column(df, ["sc1_vel_3", "sc1_velocity_3"], "velocity-z")
    t = _to_float_vector(df.time, "simulation-time")

    n = length(t)
    n >= 3 || throw(ArgumentError("Not enough rows to extract extrema (need >=3, got $n)."))

    r = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
    speed_kmps = sqrt.(vx .^ 2 .+ vy .^ 2 .+ vz .^ 2) .* 1e-3
    alt_km = Vector{Float64}(undef, n)
    if altitude_mode == :vacuum
        @inbounds for i in eachindex(r)
            alt_km[i] = (r[i] - planet.Rp_e) * 1e-3
        end
    elseif altitude_mode == :oblate
        @inbounds for i in eachindex(r)
            pos_i = SVector{3, Float64}(x[i], y[i], z[i])
            vel_i = SVector{3, Float64}(vx[i], vy[i], vz[i])
            pos_p, _ = r_intor_p!(pos_i, vel_i, planet)
            alt_km[i] = rtolatlong(pos_p, planet)[1] * 1e-3
        end
    else
        throw(ArgumentError("Unsupported altitude_mode='$altitude_mode' in _extract_extrema_series"))
    end
    rdot = (x .* vx .+ y .* vy .+ z .* vz) ./ r

    peri_t = Float64[]
    peri_alt = Float64[]
    peri_speed = Float64[]
    apo_t = Float64[]
    apo_alt = Float64[]
    apo_speed = Float64[]
    last_peri_t = -Inf
    last_apo_t = -Inf
    min_sep_s = 500.0

    for i in 1:(n - 1)
        s0 = rdot[i]
        s1 = rdot[i + 1]
        if !(isfinite(s0) && isfinite(s1))
            continue
        end

        denom = abs(s0) + abs(s1)
        w = denom == 0.0 ? 0.5 : abs(s0) / denom
        t_evt = (1.0 - w) * t[i] + w * t[i + 1]
        alt_evt = (1.0 - w) * alt_km[i] + w * alt_km[i + 1]
        speed_evt = (1.0 - w) * speed_kmps[i] + w * speed_kmps[i + 1]

        if s0 < 0.0 && s1 >= 0.0
            if isempty(peri_t) || (t_evt - last_peri_t) >= min_sep_s
                push!(peri_t, t_evt)
                push!(peri_alt, alt_evt)
                push!(peri_speed, speed_evt)
                last_peri_t = t_evt
            elseif alt_evt < peri_alt[end]
                peri_t[end] = t_evt
                peri_alt[end] = alt_evt
                peri_speed[end] = speed_evt
                last_peri_t = t_evt
            end
        elseif s0 > 0.0 && s1 <= 0.0
            if isempty(apo_t) || (t_evt - last_apo_t) >= min_sep_s
                push!(apo_t, t_evt)
                push!(apo_alt, alt_evt)
                push!(apo_speed, speed_evt)
                last_apo_t = t_evt
            elseif alt_evt > apo_alt[end]
                apo_t[end] = t_evt
                apo_alt[end] = alt_evt
                apo_speed[end] = speed_evt
                last_apo_t = t_evt
            end
        end
    end

    if isempty(peri_alt) || isempty(apo_alt)
        # Fallback for coarse/edge trajectories where radial-velocity sign changes
        # are sparse in sampled output; detect local extrema directly in altitude.
        last_peri_t = -Inf
        last_apo_t = -Inf
        for i in 2:(n - 1)
            a0 = alt_km[i - 1]
            a1 = alt_km[i]
            a2 = alt_km[i + 1]
            ti = t[i]
            vi = speed_kmps[i]
            if a1 <= a0 && a1 < a2
                if isempty(peri_t) || (ti - last_peri_t) >= min_sep_s
                    push!(peri_t, ti)
                    push!(peri_alt, a1)
                    push!(peri_speed, vi)
                    last_peri_t = ti
                elseif a1 < peri_alt[end]
                    peri_t[end] = ti
                    peri_alt[end] = a1
                    peri_speed[end] = vi
                    last_peri_t = ti
                end
            elseif a1 >= a0 && a1 > a2
                if isempty(apo_t) || (ti - last_apo_t) >= min_sep_s
                    push!(apo_t, ti)
                    push!(apo_alt, a1)
                    push!(apo_speed, vi)
                    last_apo_t = ti
                elseif a1 > apo_alt[end]
                    apo_t[end] = ti
                    apo_alt[end] = a1
                    apo_speed[end] = vi
                    last_apo_t = ti
                end
            end
        end
    end

    peri_orbit = collect(1.0:length(peri_alt))
    apo_orbit = collect(1.0:length(apo_alt))
    return (
        peri=(orbit=peri_orbit, altitude=peri_alt, time_s=peri_t, speed_kmps=peri_speed),
        apo=(orbit=apo_orbit, altitude=apo_alt, time_s=apo_t, speed_kmps=apo_speed)
    )
end

function _load_telemetry_curve(path::String, max_points::Int)
    df = DataFrame(Arrow.Table(path))
    @assert all(name in names(df) for name in ["orbit", "altitude"])
    sort!(df, :orbit)
    if max_points > 0 && nrow(df) > max_points
        df = first(df, max_points)
    end
    return (
        orbit=_to_float_vector(df.orbit, "telemetry-orbit"),
        altitude=_to_float_vector(df.altitude, "telemetry-altitude")
    )
end

function _load_time_aligned_telemetry(cfg::TimeAlignedScenarioConfig, max_points::Int)
    df = DataFrame(Arrow.Table(cfg.telemetry_path))

    time_s = _require_column(df, [cfg.telemetry_time_col], "telemetry-time")
    altitude_km = _require_column(df, [cfg.telemetry_altitude_col], "telemetry-altitude")
    x_km = _require_column(df, [cfg.telemetry_x_col], "telemetry-x")
    y_km = _require_column(df, [cfg.telemetry_y_col], "telemetry-y")
    z_km = _require_column(df, [cfg.telemetry_z_col], "telemetry-z")
    sma_km = _require_column(df, [cfg.telemetry_sma_col], "telemetry-sma")
    ecc = _require_column(df, [cfg.telemetry_ecc_col], "telemetry-ecc")
    inc_deg = _require_column(df, [cfg.telemetry_inc_col], "telemetry-inc")
    aop_deg = _require_column(df, [cfg.telemetry_aop_col], "telemetry-aop")
    raan_deg = _require_column(df, [cfg.telemetry_raan_col], "telemetry-raan")
    ta_deg = _require_column(df, [cfg.telemetry_ta_col], "telemetry-ta")

    perm = sortperm(time_s)
    time_s = time_s[perm]
    altitude_km = altitude_km[perm]
    x_km = x_km[perm]
    y_km = y_km[perm]
    z_km = z_km[perm]
    sma_km = sma_km[perm]
    ecc = ecc[perm]
    inc_deg = inc_deg[perm]
    aop_deg = aop_deg[perm]
    raan_deg = raan_deg[perm]
    ta_deg = ta_deg[perm]

    if max_points > 0 && length(time_s) > max_points
        keep = 1:max_points
        time_s = time_s[keep]
        altitude_km = altitude_km[keep]
        x_km = x_km[keep]
        y_km = y_km[keep]
        z_km = z_km[keep]
        sma_km = sma_km[keep]
        ecc = ecc[keep]
        inc_deg = inc_deg[keep]
        aop_deg = aop_deg[keep]
        raan_deg = raan_deg[keep]
        ta_deg = ta_deg[keep]
    end

    length(time_s) >= 2 || throw(ArgumentError("Need at least 2 telemetry samples for $(cfg.name), got $(length(time_s))."))
    t0 = time_s[1]
    time_s = time_s .- t0

    return (
        time_s=time_s,
        altitude_km=altitude_km,
        x_km=x_km,
        y_km=y_km,
        z_km=z_km,
        vx_kmps=_differentiate_series(x_km, time_s),
        vy_kmps=_differentiate_series(y_km, time_s),
        vz_kmps=_differentiate_series(z_km, time_s),
        sma_km=sma_km[1],
        ecc=ecc[1],
        inc_deg=inc_deg[1],
        aop_deg=aop_deg[1],
        raan_deg=raan_deg[1],
        ta_deg=ta_deg[1]
    )
end

function _differentiate_series(values::Vector{Float64}, time_s::Vector{Float64})::Vector{Float64}
    n = length(values)
    n == length(time_s) || throw(ArgumentError("values/time length mismatch in _differentiate_series: $n vs $(length(time_s))"))
    n >= 2 || throw(ArgumentError("Need at least 2 samples in _differentiate_series (got $n)."))
    dv = Vector{Float64}(undef, n)
    @inbounds begin
        dt0 = max(time_s[2] - time_s[1], eps(Float64))
        dv[1] = (values[2] - values[1]) / dt0
        for i in 2:(n - 1)
            dt = max(time_s[i + 1] - time_s[i - 1], eps(Float64))
            dv[i] = (values[i + 1] - values[i - 1]) / dt
        end
        dtn = max(time_s[n] - time_s[n - 1], eps(Float64))
        dv[n] = (values[n] - values[n - 1]) / dtn
    end
    return dv
end

function _extract_extrema_from_time_aligned_telemetry(
    time_s::Vector{Float64},
    altitude_km::Vector{Float64},
    min_sep_s::Float64;
    speed_kmps::Union{Nothing, Vector{Float64}}=nothing
)
    n = length(time_s)
    n == length(altitude_km) || throw(ArgumentError("time/altitude length mismatch: $n vs $(length(altitude_km))"))
    n >= 3 || throw(ArgumentError("Need at least 3 telemetry samples to derive peri/apo extrema (got $n)."))
    speed_series = if speed_kmps === nothing
        fill(NaN, n)
    else
        length(speed_kmps) == n || throw(ArgumentError(
            "time/speed length mismatch in _extract_extrema_from_time_aligned_telemetry: $n vs $(length(speed_kmps))"
        ))
        speed_kmps
    end

    peri_t = Float64[]
    peri_alt = Float64[]
    peri_speed = Float64[]
    apo_t = Float64[]
    apo_alt = Float64[]
    apo_speed = Float64[]
    last_peri_t = -Inf
    last_apo_t = -Inf

    @inbounds for i in 2:(n - 1)
        a0 = altitude_km[i - 1]
        a1 = altitude_km[i]
        a2 = altitude_km[i + 1]
        ti = time_s[i]
        vi = speed_series[i]
        if a1 <= a0 && a1 < a2
            if isempty(peri_t) || (ti - last_peri_t) >= min_sep_s
                push!(peri_t, ti)
                push!(peri_alt, a1)
                push!(peri_speed, vi)
                last_peri_t = ti
            elseif a1 < peri_alt[end]
                peri_t[end] = ti
                peri_alt[end] = a1
                peri_speed[end] = vi
                last_peri_t = ti
            end
        elseif a1 >= a0 && a1 > a2
            if isempty(apo_t) || (ti - last_apo_t) >= min_sep_s
                push!(apo_t, ti)
                push!(apo_alt, a1)
                push!(apo_speed, vi)
                last_apo_t = ti
            elseif a1 > apo_alt[end]
                apo_t[end] = ti
                apo_alt[end] = a1
                apo_speed[end] = vi
                last_apo_t = ti
            end
        end
    end

    isempty(peri_alt) && throw(ArgumentError("No periapsis extrema derived from telemetry altitude history."))
    isempty(apo_alt) && throw(ArgumentError("No apoapsis extrema derived from telemetry altitude history."))

    return (
        peri=(orbit=collect(1.0:length(peri_alt)), altitude=peri_alt, time_s=peri_t, speed_kmps=peri_speed),
        apo=(orbit=collect(1.0:length(apo_alt)), altitude=apo_alt, time_s=apo_t, speed_kmps=apo_speed)
    )
end

function _interp_linear(x::Vector{Float64}, y::Vector{Float64}, xq::Vector{Float64})
    length(x) == length(y) || throw(ArgumentError("x/y length mismatch: $(length(x)) vs $(length(y))"))
    n = length(x)
    n >= 1 || throw(ArgumentError("Interpolation domain is empty."))
    n >= 2 || return fill(y[1], length(xq))

    out = Vector{Float64}(undef, length(xq))
    j = 1
    @inbounds for i in eachindex(xq)
        q = xq[i]
        q <= x[1] && (out[i] = y[1]; continue)
        q >= x[end] && (out[i] = y[end]; continue)
        while j < n - 1 && x[j + 1] < q
            j += 1
        end
        x0, x1 = x[j], x[j + 1]
        y0, y1 = y[j], y[j + 1]
        w = x1 == x0 ? 0.0 : (q - x0) / (x1 - x0)
        out[i] = y0 + w * (y1 - y0)
    end
    return out
end

@inline function _normalized_axis(n::Int)
    n <= 1 && return [0.0]
    return collect(range(0.0, 1.0, length=n))
end

function _compare_orbit_curve(
    scenario::String,
    event::String,
    telemetry_axis::Vector{Float64},
    telemetry_values::Vector{Float64},
    sim_values::Vector{Float64};
    sim_axis::Union{Nothing, Vector{Float64}}=nothing
)
    n_tel = length(telemetry_values)
    n_sim = length(sim_values)
    n_tel == 0 && error("Telemetry series for $scenario/$event is empty.")

    if n_sim == 0
        return (
            scenario=scenario,
            event=event,
            n_telemetry=n_tel,
            n_sim=0,
            coverage=0.0,
            mae_km=Inf,
            rmse_km=Inf,
            max_abs_km=Inf,
            p95_abs_km=Inf,
            bias_km=Inf,
            nmae=Inf,
            nrmse=Inf,
            telemetry_axis_start=telemetry_axis[1],
            telemetry_axis_end=telemetry_axis[end]
        ), DataFrame(
            scenario=String[],
            event=String[],
            idx=Int[],
            telemetry_axis=Float64[],
            telemetry_value_km=Float64[],
            sim_interp_value_km=Float64[],
            error_km=Float64[]
        )
    end

    sim_interp = if sim_axis === nothing
        u_tel = _normalized_axis(n_tel)
        u_sim = _normalized_axis(n_sim)
        _interp_linear(u_sim, sim_values, u_tel)
    else
        length(sim_axis) == n_sim || throw(ArgumentError(
            "sim_axis length ($(length(sim_axis))) must match sim_values length ($n_sim) for $scenario/$event"
        ))
        _interp_linear(sim_axis, sim_values, telemetry_axis)
    end

    err = sim_interp .- telemetry_values
    abs_err = abs.(err)
    tel_range = max(maximum(telemetry_values) - minimum(telemetry_values), 1e-9)

    return (
        scenario=scenario,
        event=event,
        n_telemetry=n_tel,
        n_sim=n_sim,
        coverage=min(n_tel, n_sim) / n_tel,
        mae_km=mean(abs_err),
        rmse_km=sqrt(mean(err .^ 2)),
        max_abs_km=maximum(abs_err),
        p95_abs_km=quantile(abs_err, 0.95),
        bias_km=mean(err),
        nmae=mean(abs_err) / tel_range,
        nrmse=sqrt(mean(err .^ 2)) / tel_range,
        telemetry_axis_start=telemetry_axis[1],
        telemetry_axis_end=telemetry_axis[end]
    ), DataFrame(
        scenario=fill(scenario, n_tel),
        event=fill(event, n_tel),
        idx=collect(1:n_tel),
        telemetry_axis=telemetry_axis,
        telemetry_value_km=telemetry_values,
        sim_interp_value_km=sim_interp,
        error_km=err
    )
end

function _compare_time_series(
    scenario::String,
    event::String,
    telemetry_time::Vector{Float64},
    telemetry_values::Vector{Float64},
    sim_time::Vector{Float64},
    sim_values::Vector{Float64}
)
    n_tel = length(telemetry_values)
    n_sim = length(sim_values)
    n_tel == 0 && error("Telemetry series for $scenario/$event is empty.")

    if n_sim == 0
        return (
            scenario=scenario,
            event=event,
            n_telemetry=n_tel,
            n_sim=0,
            coverage=0.0,
            mae_km=Inf,
            rmse_km=Inf,
            max_abs_km=Inf,
            p95_abs_km=Inf,
            bias_km=Inf,
            nmae=Inf,
            nrmse=Inf,
            telemetry_axis_start=telemetry_time[1],
            telemetry_axis_end=telemetry_time[end]
        ), DataFrame(
            scenario=String[],
            event=String[],
            idx=Int[],
            telemetry_axis=Float64[],
            telemetry_value_km=Float64[],
            sim_interp_value_km=Float64[],
            error_km=Float64[]
        )
    end

    sim_interp = _interp_linear(sim_time, sim_values, telemetry_time)
    err = sim_interp .- telemetry_values
    abs_err = abs.(err)
    tel_range = max(maximum(telemetry_values) - minimum(telemetry_values), 1e-9)

    return (
        scenario=scenario,
        event=event,
        n_telemetry=n_tel,
        n_sim=n_sim,
        coverage=min(n_tel, n_sim) / n_tel,
        mae_km=mean(abs_err),
        rmse_km=sqrt(mean(err .^ 2)),
        max_abs_km=maximum(abs_err),
        p95_abs_km=quantile(abs_err, 0.95),
        bias_km=mean(err),
        nmae=mean(abs_err) / tel_range,
        nrmse=sqrt(mean(err .^ 2)) / tel_range,
        telemetry_axis_start=telemetry_time[1],
        telemetry_axis_end=telemetry_time[end]
    ), DataFrame(
        scenario=fill(scenario, n_tel),
        event=fill(event, n_tel),
        idx=collect(1:n_tel),
        telemetry_axis=telemetry_time,
        telemetry_value_km=telemetry_values,
        sim_interp_value_km=sim_interp,
        error_km=err
    )
end

@inline function _calibration_active(cal::CalibrationConfig, profile::Symbol)::Bool
    return cal.enabled && (profile in cal.profiles)
end

@inline function _grid_values(min_v::Float64, max_v::Float64, steps::Int)::Vector{Float64}
    if steps <= 1 || isapprox(min_v, max_v; rtol=0.0, atol=1e-12)
        return [min_v]
    end
    return collect(range(min_v, max_v, length=steps))
end

function _estimate_event_biases(error_tables::Vector{DataFrame}, bias_abs_max_km::Float64)::Dict{String, Float64}
    out = Dict{String, Float64}()
    for df in error_tables
        nrow(df) == 0 && continue
        event = String(df.event[1])
        bias = -mean(_to_float_vector(df.error_km, "error_km:$event"))
        out[event] = clamp(bias, -bias_abs_max_km, bias_abs_max_km)
    end
    return out
end

function _calibration_score(rows::AbstractVector{<:NamedTuple}, objective::String)::Float64
    isempty(rows) && return Inf
    if objective == "mean_nmae"
        return mean([Float64(r.nmae) for r in rows])
    elseif objective == "mean_rmse_km"
        return mean([Float64(r.rmse_km) for r in rows])
    elseif objective == "max_nmae"
        return maximum([Float64(r.nmae) for r in rows])
    end
    throw(ArgumentError("Unsupported calibration objective '$objective'"))
end

function _annotate_calibration_rows(
    rows::AbstractVector{<:NamedTuple},
    cd_scale::Float64,
    cr_value::Float64,
    bias_by_event::Dict{String, Float64},
    score::Float64,
    selected_runtime_s::Float64,
    dt_max_orbit_s::Float64,
    calibration_runtime_s::Float64,
    calibration_used::Bool,
    solver_info
)::Vector{NamedTuple}
    out = NamedTuple[]
    for row in rows
        push!(
            out,
            merge(
                row,
                (
                    calibration_used=calibration_used,
                    calibrated_cd_scale=cd_scale,
                    calibrated_cr=cr_value,
                    calibrated_bias_km=get(bias_by_event, String(row.event), 0.0),
                    calibration_score=score,
                    selected_simulation_runtime_s=selected_runtime_s,
                    dt_max_orbit_s=dt_max_orbit_s,
                    calibration_runtime_s=calibration_runtime_s,
                    solver_mode=solver_info.solver_mode,
                    solver_sequence=solver_info.solver_sequence,
                    solver_fallback_used=solver_info.solver_fallback_used,
                    solver_fallback_count=solver_info.solver_fallback_count,
                    solver_fallback_trigger=solver_info.solver_fallback_trigger,
                    solver_retcode=solver_info.solver_retcode,
                    solver_maxiters=solver_info.solver_maxiters,
                    solver_maxiters_retry_used=solver_info.solver_maxiters_retry_used
                )
            )
        )
    end
    return out
end

function _orbit_rows_errors(
    cfg::OrbitEventsScenarioConfig,
    args::SimulationConfiguration,
    results_df::DataFrame,
    max_points::Int;
    bias_by_event::Dict{String, Float64}=Dict{String, Float64}()
)
    extrema = _extract_extrema_series(results_df, args.environment_model.planet, cfg.orbit_altitude_mode)
    tele_peri = _load_telemetry_curve(cfg.telemetry_peri_path, max_points)
    tele_apo = _load_telemetry_curve(cfg.telemetry_apo_path, max_points)
    peri_bias = get(bias_by_event, "peri", 0.0)
    apo_bias = get(bias_by_event, "apo", 0.0)
    peri_step = length(tele_peri.orbit) >= 2 ? median(diff(tele_peri.orbit)) : 1.0
    apo_step = length(tele_apo.orbit) >= 2 ? median(diff(tele_apo.orbit)) : 1.0
    peri_sim_axis = tele_peri.orbit[1] .+ peri_step .* collect(0:(length(extrema.peri.altitude)-1))
    apo_sim_axis = tele_apo.orbit[1] .+ apo_step .* collect(0:(length(extrema.apo.altitude)-1))
    peri_summary, peri_errors = _compare_orbit_curve(
        cfg.name,
        "peri",
        tele_peri.orbit,
        tele_peri.altitude,
        extrema.peri.altitude .+ peri_bias;
        sim_axis=peri_sim_axis
    )
    apo_summary, apo_errors = _compare_orbit_curve(
        cfg.name,
        "apo",
        tele_apo.orbit,
        tele_apo.altitude,
        extrema.apo.altitude .+ apo_bias;
        sim_axis=apo_sim_axis
    )
    return [peri_summary, apo_summary], [peri_errors, apo_errors]
end

function _time_aligned_rows_errors(
    cfg::TimeAlignedScenarioConfig,
    args::SimulationConfiguration,
    results_df::DataFrame,
    telemetry;
    bias_by_event::Dict{String, Float64}=Dict{String, Float64}()
)
    if cfg.comparison_mode == :orbit_events
        tele_speed_kmps = sqrt.(telemetry.vx_kmps .^ 2 .+ telemetry.vy_kmps .^ 2 .+ telemetry.vz_kmps .^ 2)
        tele_extrema = _extract_extrema_from_time_aligned_telemetry(
            telemetry.time_s,
            telemetry.altitude_km,
            cfg.extrema_min_separation_s;
            speed_kmps=tele_speed_kmps
        )
        sim_extrema = _extract_extrema_series(results_df, args.environment_model.planet, cfg.orbit_altitude_mode)
        peri_bias = get(bias_by_event, "peri", 0.0)
        apo_bias = get(bias_by_event, "apo", 0.0)
        peri_speed_bias = get(bias_by_event, "peri_speed", 0.0)
        apo_speed_bias = get(bias_by_event, "apo_speed", 0.0)

        peri_summary, peri_errors = _compare_orbit_curve(
            cfg.name,
            "peri",
            tele_extrema.peri.orbit,
            tele_extrema.peri.altitude,
            sim_extrema.peri.altitude .+ peri_bias
        )
        apo_summary, apo_errors = _compare_orbit_curve(
            cfg.name,
            "apo",
            tele_extrema.apo.orbit,
            tele_extrema.apo.altitude,
            sim_extrema.apo.altitude .+ apo_bias
        )

        peri_speed_summary, peri_speed_errors = _compare_orbit_curve(
            cfg.name,
            "peri_speed",
            tele_extrema.peri.orbit,
            tele_extrema.peri.speed_kmps,
            sim_extrema.peri.speed_kmps .+ peri_speed_bias
        )
        apo_speed_summary, apo_speed_errors = _compare_orbit_curve(
            cfg.name,
            "apo_speed",
            tele_extrema.apo.orbit,
            tele_extrema.apo.speed_kmps,
            sim_extrema.apo.speed_kmps .+ apo_speed_bias
        )

        return [peri_summary, apo_summary, peri_speed_summary, apo_speed_summary], [peri_errors, apo_errors, peri_speed_errors, apo_speed_errors]
    end

    sim_time = _to_float_vector(results_df.time, "sim-time")
    sim_x_m = _require_column(results_df, ["sc1_pos_1", "sc1_position_1"], "sim-position-x")
    sim_y_m = _require_column(results_df, ["sc1_pos_2", "sc1_position_2"], "sim-position-y")
    sim_z_m = _require_column(results_df, ["sc1_pos_3", "sc1_position_3"], "sim-position-z")

    sim_x_km = sim_x_m .* 1e-3 .+ get(bias_by_event, "state_x_time", 0.0)
    sim_y_km = sim_y_m .* 1e-3 .+ get(bias_by_event, "state_y_time", 0.0)
    sim_z_km = sim_z_m .* 1e-3 .+ get(bias_by_event, "state_z_time", 0.0)
    sim_r_m = sqrt.(sim_x_m .^ 2 .+ sim_y_m .^ 2 .+ sim_z_m .^ 2)
    sim_altitude_km = (sim_r_m .- args.environment_model.planet.Rp_e) .* 1e-3 .+ get(bias_by_event, "altitude_time", 0.0)

    altitude_summary, altitude_errors = _compare_time_series(
        cfg.name,
        "altitude_time",
        telemetry.time_s,
        telemetry.altitude_km,
        sim_time,
        sim_altitude_km
    )
    x_summary, x_errors = _compare_time_series(
        cfg.name,
        "state_x_time",
        telemetry.time_s,
        telemetry.x_km,
        sim_time,
        sim_x_km
    )
    y_summary, y_errors = _compare_time_series(
        cfg.name,
        "state_y_time",
        telemetry.time_s,
        telemetry.y_km,
        sim_time,
        sim_y_km
    )
    z_summary, z_errors = _compare_time_series(
        cfg.name,
        "state_z_time",
        telemetry.time_s,
        telemetry.z_km,
        sim_time,
        sim_z_km
    )
    return [altitude_summary, x_summary, y_summary, z_summary], [altitude_errors, x_errors, y_errors, z_errors]
end

@inline _tolerances_for(cfg::OrbitEventsScenarioConfig, profile::Symbol) = profile == :quick ? cfg.tolerances_quick : cfg.tolerances_full
@inline _tolerances_for(cfg::TimeAlignedScenarioConfig, profile::Symbol) = profile == :quick ? cfg.tolerances_quick : cfg.tolerances_full

@inline function _default_plots_outdir(out_summary::String, profile::Symbol)::String
    return normpath(joinpath(dirname(out_summary), "telemetry_plots_$(String(profile))"))
end

function _generate_plots(summary_csv::String, errors_csv::String, profile::Symbol)::String
    plot_script = joinpath(REPO_ROOT, "src", "analysis", "plotting", "telemetry_orbit_accuracy_plots.jl")
    isfile(plot_script) || throw(ErrorException("Missing plotting script: $plot_script"))
    outdir = _default_plots_outdir(summary_csv, profile)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(joinpath(REPO_ROOT, ".AGORA")) $plot_script --summary=$summary_csv --errors=$errors_csv --outdir=$outdir`
    run(cmd)
    return outdir
end

@inline _axis_units(cfg::OrbitEventsScenarioConfig) = cfg.units_x
@inline _axis_units(cfg::TimeAlignedScenarioConfig) = cfg.units_x

@inline function _value_units(cfg::OrbitEventsScenarioConfig, event::String)
    if occursin("speed", event) || occursin("vx", event) || occursin("vy", event) || occursin("vz", event)
        return get(cfg.units_y, event, "km/s")
    end
    return get(cfg.units_y, event, "km")
end
@inline function _value_units(cfg::TimeAlignedScenarioConfig, event::String)
    if occursin("speed", event) || occursin("vx", event) || occursin("vy", event) || occursin("vz", event)
        return get(cfg.units_y, event, "km/s")
    end
    return get(cfg.units_y, event, "km")
end

@inline function _display_value_units(value_units::AbstractString)::String
    token = lowercase(strip(String(value_units)))
    return token == "km/s" ? "m/s" : String(value_units)
end

@inline function _display_value_scale(value_units::AbstractString)::Float64
    token = lowercase(strip(String(value_units)))
    return token == "km/s" ? 1e3 : 1.0
end

function _append_display_metric_columns!(summary_df::DataFrame)::Nothing
    nrow(summary_df) == 0 && return nothing
    hasproperty(summary_df, :value_units) || return nothing

    units_raw = [String(v) for v in summary_df.value_units]
    scales = [_display_value_scale(u) for u in units_raw]
    units_display = [_display_value_units(u) for u in units_raw]
    summary_df.value_units_display = units_display

    function _add_scaled_column!(src::Symbol, dst::Symbol)
        hasproperty(summary_df, src) || return nothing
        src_values = summary_df[!, src]
        out = Vector{Float64}(undef, length(src_values))
        @inbounds for i in eachindex(src_values)
            out[i] = Float64(src_values[i]) * scales[i]
        end
        summary_df[!, dst] = out
        return nothing
    end

    _add_scaled_column!(:mae_km, :mae_display)
    _add_scaled_column!(:rmse_km, :rmse_display)
    _add_scaled_column!(:max_abs_km, :max_abs_display)
    _add_scaled_column!(:p95_abs_km, :p95_abs_display)
    _add_scaled_column!(:bias_km, :bias_display)
    _add_scaled_column!(:limit_max_abs_km, :limit_max_abs_display)
    _add_scaled_column!(:limit_max_rmse_km, :limit_max_rmse_display)
    return nothing
end

function _append_display_error_columns!(errors_df::DataFrame, summary_df::DataFrame)::Nothing
    nrow(errors_df) == 0 && return nothing
    hasproperty(errors_df, :scenario) || return nothing
    hasproperty(errors_df, :event) || return nothing

    unit_map = Dict{Tuple{String, String}, Tuple{Float64, String}}()
    if nrow(summary_df) > 0 && hasproperty(summary_df, :scenario) && hasproperty(summary_df, :event) && hasproperty(summary_df, :value_units)
        for row in eachrow(summary_df)
            key = (String(row.scenario), String(row.event))
            units = String(row.value_units)
            unit_map[key] = (_display_value_scale(units), _display_value_units(units))
        end
    end

    n = nrow(errors_df)
    scales = Vector{Float64}(undef, n)
    units_display = Vector{String}(undef, n)
    @inbounds for i in 1:n
        key = (String(errors_df.scenario[i]), String(errors_df.event[i]))
        scale, unit = get(unit_map, key, (1.0, "km"))
        scales[i] = scale
        units_display[i] = unit
    end
    errors_df.value_units_display = units_display

    function _add_scaled_column!(src::Symbol, dst::Symbol)
        hasproperty(errors_df, src) || return nothing
        src_values = errors_df[!, src]
        out = Vector{Float64}(undef, length(src_values))
        @inbounds for i in eachindex(src_values)
            out[i] = Float64(src_values[i]) * scales[i]
        end
        errors_df[!, dst] = out
        return nothing
    end

    _add_scaled_column!(:telemetry_value_km, :telemetry_value_display)
    _add_scaled_column!(:sim_interp_value_km, :sim_interp_value_display)
    _add_scaled_column!(:error_km, :error_display)
    return nothing
end

@inline _source_file(cfg::OrbitEventsScenarioConfig, event::String) = event == "peri" ? cfg.telemetry_peri_path : cfg.telemetry_apo_path
@inline _source_file(cfg::TimeAlignedScenarioConfig, event::String) = cfg.telemetry_path

@inline _orbit_altitude_mode(cfg::OrbitEventsScenarioConfig) = String(cfg.orbit_altitude_mode)
@inline function _orbit_altitude_mode(cfg::TimeAlignedScenarioConfig)
    return cfg.comparison_mode == :orbit_events ? String(cfg.orbit_altitude_mode) : "n/a"
end

@inline _maneuver_count(cfg::OrbitEventsScenarioConfig) = length(cfg.maneuver_orbit_numbers)
@inline _maneuver_count(::TimeAlignedScenarioConfig) = 0

@inline _min_eval_points(cfg::OrbitEventsScenarioConfig) = cfg.min_eval_points
@inline _min_eval_points(cfg::TimeAlignedScenarioConfig) = cfg.min_eval_points

@inline function _scenario_status_extra(cfg::OrbitEventsScenarioConfig)::String
    return ", altitude_mode=$(String(cfg.orbit_altitude_mode)), maneuvers=$(length(cfg.maneuver_orbit_numbers))"
end
@inline function _scenario_status_extra(cfg::TimeAlignedScenarioConfig)::String
    if cfg.comparison_mode == :orbit_events
        return ", comparison_mode=orbit_events, altitude_mode=$(String(cfg.orbit_altitude_mode))"
    end
    return ", comparison_mode=time_aligned_state"
end

function _evaluate_thresholds(row, cfg::AbstractScenarioConfig, profile::Symbol)
    tmap = _tolerances_for(cfg, profile)
    event_name = String(row.event)
    tol = get(tmap, event_name, nothing)
    if tol === nothing && endswith(event_name, "_speed")
        base_event = replace(event_name, "_speed" => "")
        tol = get(tmap, base_event, nothing)
    end
    tol === nothing && throw(ArgumentError("No tolerance configured for event '$(row.event)' in scenario '$(row.scenario)'"))

    pass_points = row.n_sim >= _min_eval_points(cfg)
    pass_abs = row.max_abs_km <= tol.max_abs_km
    pass_nmae = row.nmae <= tol.max_nmae
    pass_rmse = row.rmse_km <= tol.max_rmse_km
    pass = pass_points && pass_abs && pass_nmae && pass_rmse

    return (
        pass=pass,
        pass_points=pass_points,
        pass_abs=pass_abs,
        pass_nmae=pass_nmae,
        pass_rmse=pass_rmse,
        limit_max_abs_km=tol.max_abs_km,
        limit_nmae=tol.max_nmae,
        limit_max_rmse_km=tol.max_rmse_km,
        min_eval_points=_min_eval_points(cfg)
    )
end

function _run_simulation_dataframe(args::SimulationConfiguration, scenario_name::String, truth::AtmosphereTruthConfig, profile::Symbol)
    return mktempdir() do tmp
        cfg_run = SimulationConfiguration(
            file_paths=args.file_paths,
            simulation_settings=SimulationSettings(
                results=true,
                verbose=false,
                results_directory=tmp,
                generate_plots=false,
                generate_filenames=false,
                normalize=false,
                save_csv=true
            ),
            mission_configuration=args.mission_configuration,
            environment_model=args.environment_model,
            dynamics_model=args.dynamics_model,
            guidance_model=args.guidance_model,
            navigation_model=args.navigation_model,
            control_model=args.control_model,
            initial_time=args.initial_time,
            integration_tolerances=args.integration_tolerances
        )

        save_fields = _save_fields_for_study()
        base_maxiters = _telemetry_solver_maxiters(profile)

        function _run_once(maxiters::Int)
            solve_result = nothing
            elapsed_s = @elapsed begin
                solver_mode = _telemetry_solver_mode()
                withenv(
                    "SPACEAGORA_WARN_NORMALIZE" => "0",
                    "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
                    "SPACEAGORA_SOLVER_MODE" => solver_mode,
                    "SPACEAGORA_SOLVER_MAXITERS" => string(maxiters),
                    "SPACEAGORA_GRAM_OFFLINE_SURROGATE" => truth.gram_offline_surrogate,
                    "SPACEAGORA_GRAM_STATIC_GRID" => truth.gram_static_grid ? "on" : "off",
                    "SPACEAGORA_GRAM_TRACK_CACHE" => truth.gram_track_cache ? "on" : "off",
                    "SPACEAGORA_GRAM_GLOBAL_LOCK" => truth.gram_global_lock
                ) do
                    cd(tmp) do
                        solve_result = SimulationEngine.run_simulation(
                            cfg_run;
                            isolate_state=false,
                            save_fields=save_fields,
                            return_solution=true,
                            return_solver_metadata=true
                        )
                    end
                end
            end
            return solve_result, elapsed_s
        end

        retry_used = false
        maxiters_used = base_maxiters
        solve_result = nothing
        elapsed_s = 0.0
        try
            solve_result, elapsed_s = _run_once(base_maxiters)
        catch err
            if _is_maxiters_error(err)
                retry_used = true
                maxiters_used = _telemetry_solver_retry_maxiters(base_maxiters)
                @warn "Telemetry scenario '$scenario_name' hit MaxIters; retrying with larger SPACEAGORA_SOLVER_MAXITERS." base_maxiters=base_maxiters retry_maxiters=maxiters_used
                solve_result, elapsed_s = _run_once(maxiters_used)
            else
                rethrow(err)
            end
        end

        csv_path = joinpath(tmp, "simulation_results.csv")
        isfile(csv_path) || error("Missing simulation output for $scenario_name: $csv_path")
        results_df = CSV.read(csv_path, DataFrame)

        solver_trace = solve_result.solver_trace
        solver_sequence = isempty(solver_trace) ? "n/a" : join([meta.solver for meta in solver_trace], "->")
        solver_fallback_count = count(meta -> meta.fallback_used, solver_trace)
        fallback_triggers = [meta.trigger_retcode for meta in solver_trace if !(meta.trigger_retcode isa Missing)]
        solver_fallback_trigger = isempty(fallback_triggers) ? "n/a" : join(sort(unique(string.(fallback_triggers))), "|")
        solver_info = (
            solver_mode=String(solve_result.solver_mode),
            solver_sequence=solver_sequence,
            solver_fallback_used=solver_fallback_count > 0,
            solver_fallback_count=solver_fallback_count,
            solver_fallback_trigger=solver_fallback_trigger,
            solver_retcode=string(solve_result.solution.retcode),
            solver_maxiters=maxiters_used,
            solver_maxiters_retry_used=retry_used
        )
        return (results_df=results_df, elapsed_s=elapsed_s, solver_info=solver_info)
    end
end

function _initial_condition_from_time_aligned_telemetry(telemetry)
    sma_m = telemetry.sma_km * 1e3
    ra_m = sma_m * (1.0 + telemetry.ecc)
    rp_m = sma_m * (1.0 - telemetry.ecc)
    return InitialCondition(
        ra=ra_m,
        rp=rp_m,
        i=telemetry.inc_deg,
        ω=telemetry.aop_deg,
        Ω=telemetry.raan_deg,
        ν=telemetry.ta_deg
    )
end

function _run_single_scenario(cfg::OrbitEventsScenarioConfig, profile::Symbol)
    scenario_start = time()
    cal = cfg.calibration
    use_calibration = _calibration_active(cal, profile)

    final_is_quick = profile == :quick
    final_orbits = final_is_quick ? cfg.target_orbits_quick : cfg.target_orbits_full
    final_points = final_is_quick ? cfg.compare_points_quick : cfg.compare_points_full

    eval_profile = (use_calibration && cal.search_on_quick_subset && profile == :full) ? :quick : profile
    eval_is_quick = eval_profile == :quick
    eval_orbits = eval_is_quick ? cfg.target_orbits_quick : cfg.target_orbits_full
    eval_points = eval_is_quick ? cfg.compare_points_quick : cfg.compare_points_full

    cd_candidates = (use_calibration && cal.fit_cd_scale) ?
        _grid_values(cal.cd_scale_min, cal.cd_scale_max, cal.cd_scale_steps) : [1.0]
    cr_base = cfg.srp_cr
    cr_candidates = (use_calibration && cal.fit_cr && cfg.srp_enabled) ?
        _grid_values(cal.cr_min, cal.cr_max, cal.cr_steps) : [cr_base]

    best_cd = cd_candidates[1]
    best_cr = cr_candidates[1]
    best_score = Inf

    if use_calibration
        for cd_scale in cd_candidates, cr_value in cr_candidates
            args_eval = _make_orbit_args(cfg, eval_orbits; cd_scale=cd_scale, cr_override=cr_value)
            args_eval = _with_study_settings(args_eval; quick=eval_is_quick)
            eval_run = _run_simulation_dataframe(args_eval, cfg.name, cfg.atmosphere_truth, eval_profile)
            eval_df = eval_run.results_df
            eval_rows, eval_errors = _orbit_rows_errors(cfg, args_eval, eval_df, eval_points)
            if cal.fit_bias
                eval_bias = _estimate_event_biases(eval_errors, cal.bias_abs_max_km)
                eval_rows, eval_errors = _orbit_rows_errors(cfg, args_eval, eval_df, eval_points; bias_by_event=eval_bias)
            end
            score = _calibration_score(eval_rows, cal.objective)
            if score < best_score
                best_score = score
                best_cd = cd_scale
                best_cr = cr_value
            end
        end
    else
        best_score = NaN
    end

    args_final = _make_orbit_args(cfg, final_orbits; cd_scale=best_cd, cr_override=best_cr)
    args_final = _with_study_settings(args_final; quick=final_is_quick)
    final_run = _run_simulation_dataframe(args_final, cfg.name, cfg.atmosphere_truth, profile)
    final_df = final_run.results_df
    selected_runtime_s = final_run.elapsed_s
    solver_info = final_run.solver_info
    final_rows, final_errors = _orbit_rows_errors(cfg, args_final, final_df, final_points)
    bias_by_event = (use_calibration && cal.fit_bias) ? _estimate_event_biases(final_errors, cal.bias_abs_max_km) : Dict{String, Float64}()
    if !isempty(bias_by_event)
        final_rows, final_errors = _orbit_rows_errors(cfg, args_final, final_df, final_points; bias_by_event=bias_by_event)
    end
    final_score = use_calibration ? _calibration_score(final_rows, cal.objective) : NaN
    calibration_runtime_s = time() - scenario_start
    annotated_rows = _annotate_calibration_rows(
        final_rows,
        best_cd,
        best_cr,
        bias_by_event,
        final_score,
        selected_runtime_s,
        args_final.integration_tolerances.dt_max_orbit,
        calibration_runtime_s,
        use_calibration,
        solver_info
    )
    return annotated_rows, final_errors, calibration_runtime_s
end

function _run_single_scenario(cfg::TimeAlignedScenarioConfig, profile::Symbol)
    scenario_start = time()
    cal = cfg.calibration
    use_calibration = _calibration_active(cal, profile)

    final_is_quick = profile == :quick
    final_points = final_is_quick ? cfg.max_points_quick : cfg.max_points_full
    final_telemetry = _load_time_aligned_telemetry(cfg, final_points)
    final_ic = _initial_condition_from_time_aligned_telemetry(final_telemetry)
    final_mission_time_s = max(final_telemetry.time_s[end] - final_telemetry.time_s[1], 1.0)

    eval_profile = (use_calibration && cal.search_on_quick_subset && profile == :full) ? :quick : profile
    eval_is_quick = eval_profile == :quick
    eval_points = eval_is_quick ? cfg.max_points_quick : cfg.max_points_full
    eval_telemetry = _load_time_aligned_telemetry(cfg, eval_points)
    eval_ic = _initial_condition_from_time_aligned_telemetry(eval_telemetry)
    eval_mission_time_s = max(eval_telemetry.time_s[end] - eval_telemetry.time_s[1], 1.0)

    cd_candidates = (use_calibration && cal.fit_cd_scale) ?
        _grid_values(cal.cd_scale_min, cal.cd_scale_max, cal.cd_scale_steps) : [1.0]
    cr_base = cfg.srp_cr
    cr_candidates = (use_calibration && cal.fit_cr && cfg.srp_enabled) ?
        _grid_values(cal.cr_min, cal.cr_max, cal.cr_steps) : [cr_base]

    best_cd = cd_candidates[1]
    best_cr = cr_candidates[1]
    best_score = Inf

    if use_calibration
        for cd_scale in cd_candidates, cr_value in cr_candidates
            args_eval = _make_time_aligned_args(
                cfg,
                eval_mission_time_s,
                eval_ic;
                cd_scale=cd_scale,
                cr_override=cr_value
            )
            args_eval = _with_study_settings(args_eval; quick=eval_is_quick)
            eval_run = _run_simulation_dataframe(args_eval, cfg.name, cfg.atmosphere_truth, eval_profile)
            eval_df = eval_run.results_df
            eval_rows, eval_errors = _time_aligned_rows_errors(cfg, args_eval, eval_df, eval_telemetry)
            if cal.fit_bias
                eval_bias = _estimate_event_biases(eval_errors, cal.bias_abs_max_km)
                eval_rows, eval_errors = _time_aligned_rows_errors(cfg, args_eval, eval_df, eval_telemetry; bias_by_event=eval_bias)
            end
            score = _calibration_score(eval_rows, cal.objective)
            if score < best_score
                best_score = score
                best_cd = cd_scale
                best_cr = cr_value
            end
        end
    else
        best_score = NaN
    end

    args_final = _make_time_aligned_args(
        cfg,
        final_mission_time_s,
        final_ic;
        cd_scale=best_cd,
        cr_override=best_cr
    )
    args_final = _with_study_settings(args_final; quick=final_is_quick)
    final_run = _run_simulation_dataframe(args_final, cfg.name, cfg.atmosphere_truth, profile)
    final_df = final_run.results_df
    selected_runtime_s = final_run.elapsed_s
    solver_info = final_run.solver_info
    final_rows, final_errors = _time_aligned_rows_errors(cfg, args_final, final_df, final_telemetry)
    bias_by_event = (use_calibration && cal.fit_bias) ? _estimate_event_biases(final_errors, cal.bias_abs_max_km) : Dict{String, Float64}()
    if !isempty(bias_by_event)
        final_rows, final_errors = _time_aligned_rows_errors(cfg, args_final, final_df, final_telemetry; bias_by_event=bias_by_event)
    end
    final_score = use_calibration ? _calibration_score(final_rows, cal.objective) : NaN
    calibration_runtime_s = time() - scenario_start
    annotated_rows = _annotate_calibration_rows(
        final_rows,
        best_cd,
        best_cr,
        bias_by_event,
        final_score,
        selected_runtime_s,
        args_final.integration_tolerances.dt_max_orbit,
        calibration_runtime_s,
        use_calibration,
        solver_info
    )
    return annotated_rows, final_errors, calibration_runtime_s
end

function _run_verification(cfg::StudyConfig)::VerificationResult
    scenarios = _load_scenarios_from_manifest(cfg.manifest_path)

    summary_rows = NamedTuple[]
    error_tables = DataFrame[]
    wall_start = time()

    println("Telemetry Orbit Accuracy Study")
    println(@sprintf("profile=%s enforce=%s", String(cfg.profile), string(cfg.enforce)))
    println("manifest=$(cfg.manifest_path)")
    println("deterministic_mode=GRAM(per-scenario atmosphere_truth from manifest)")

    for sc in scenarios
        truth = sc.atmosphere_truth
        cal = sc.calibration
        println(
            "Running scenario $(sc.name) ... " *
            "(truth=$(truth.assumption_id), atmosphere=$(truth.atmosphere_model), " *
            "dataset=$(truth.atmosphere_dataset), weather=$(truth.space_weather_model), solar=$(truth.solar_flux_model), " *
            "seed=$(truth.gram_seed), " *
            "surrogate=$(truth.gram_offline_surrogate), static_grid=$(truth.gram_static_grid), " *
            "track_cache=$(truth.gram_track_cache), calibration=$(_calibration_active(cal, cfg.profile))$(_scenario_status_extra(sc)))"
        )
        row_batch, err_batch, elapsed_s = _run_single_scenario(sc, cfg.profile)
        println(@sprintf("COMPUTATIONAL TIME [%s/%s] = %.3f s", sc.name, String(cfg.profile), elapsed_s))

        for row in row_batch
            gates = _evaluate_thresholds(row, sc, cfg.profile)
            push!(summary_rows, merge(
                row,
                gates,
                (
                    profile=String(cfg.profile),
                    source_file=_source_file(sc, row.event),
                    axis_units=_axis_units(sc),
                    value_units=_value_units(sc, row.event),
                    orbit_altitude_mode=_orbit_altitude_mode(sc),
                    maneuver_count=_maneuver_count(sc),
                    simulation_runtime_s=elapsed_s,
                    timestamp_utc=string(now(UTC)),
                    atmosphere_truth_id=sc.atmosphere_truth.assumption_id,
                    atmosphere_model=sc.atmosphere_truth.atmosphere_model,
                    atmosphere_dataset=sc.atmosphere_truth.atmosphere_dataset,
                    space_weather_model=sc.atmosphere_truth.space_weather_model,
                    solar_flux_model=sc.atmosphere_truth.solar_flux_model,
                    gram_seed=sc.atmosphere_truth.gram_seed,
                    gram_perturbation_scales=join(sc.atmosphere_truth.gram_perturbation_scales, ","),
                    gram_offline_surrogate=sc.atmosphere_truth.gram_offline_surrogate,
                    gram_static_grid=sc.atmosphere_truth.gram_static_grid,
                    gram_track_cache=sc.atmosphere_truth.gram_track_cache,
                    gram_global_lock=sc.atmosphere_truth.gram_global_lock
                )
            ))
        end
        dt_by_event = Dict{String, Float64}()
        for row in row_batch
            dt_by_event[String(row.event)] = Float64(row.dt_max_orbit_s)
        end
        for err_df in err_batch
            nrow(err_df) == 0 && continue
            evt = String(err_df.event[1])
            err_df.dt_max_orbit_s = fill(get(dt_by_event, evt, NaN), nrow(err_df))
        end
        append!(error_tables, err_batch)
    end

    total_runtime_s = time() - wall_start
    println(@sprintf("TOTAL COMPUTATIONAL TIME [%s] = %.3f s", String(cfg.profile), total_runtime_s))

    summary_df = DataFrame(summary_rows)
    errors_df = isempty(error_tables) ? DataFrame() : vcat(error_tables...; cols=:union)
    if nrow(summary_df) > 0
        sort!(summary_df, [:scenario, :event])
        summary_df.total_runtime_s = fill(total_runtime_s, nrow(summary_df))
    end
    _append_display_metric_columns!(summary_df)
    _append_display_error_columns!(errors_df, summary_df)

    mkpath(dirname(cfg.out_summary))
    mkpath(dirname(cfg.out_errors))
    CSV.write(cfg.out_summary, summary_df)
    CSV.write(cfg.out_errors, errors_df)

    println("\nSummary:")
    show(summary_df, allrows=true, allcols=true)
    println("\n\nSaved:")
    println("  summary: $(cfg.out_summary)")
    println("  errors : $(cfg.out_errors)")
    plots_dir = ""
    if cfg.generate_plots
        plots_dir = _generate_plots(cfg.out_summary, cfg.out_errors, cfg.profile)
        println("  plots  : $(plots_dir)")
    else
        println("  plots  : disabled")
    end

    if cfg.enforce && nrow(summary_df) > 0
        failed = summary_df[summary_df.pass .== false, :]
        if nrow(failed) > 0
            println("\nThreshold failures:")
            show(failed, allrows=true, allcols=true)
            error("Telemetry orbit accuracy thresholds failed for $(nrow(failed)) row(s).")
        end
    end

    return VerificationResult(
        summary=summary_df,
        errors=errors_df,
        summary_path=cfg.out_summary,
        errors_path=cfg.out_errors,
        plots_dir=plots_dir,
        profile=cfg.profile,
        enforce=cfg.enforce,
        total_runtime_s=total_runtime_s
    )
end

function run_verification(request::VerificationRequest)::VerificationResult
    return _run_verification(_study_config(request))
end

function run_verification_cli(args::Vector{String}=copy(ARGS))::VerificationResult
    cfg = parse_cli(args)
    request = _request_from_study_config(cfg)
    request = VerificationRequest(
        profile=request.profile,
        out_summary=request.out_summary,
        out_errors=request.out_errors,
        manifest_path=request.manifest_path,
        enforce=cfg.enforce,
        generate_plots=cfg.generate_plots
    )
    return run_verification(request)
end

function run_study()
    result = run_verification_cli(copy(ARGS))
    return (summary=result.summary, errors=result.errors)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_verification_cli(copy(ARGS))
end

end # module TelemetryVerification
