include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl")) # Get the reference system types for the callback

using DifferentialEquations
using DiffEqBase
using LinearAlgebra
using SPICE
using Dates
using GRAMSuite
using ...RuntimeServices: SPICE_LOCK, GRAM_LOCK
using ..SimulationModel: PlanetFrameEphemerisCache, rot
using ..SimulationModel: ephemerides_time_seconds, planet_frame_lpi, ephemerides_requires_spice
using ..ParallelPolicy
using ..EnvironmentModels
using ..EnvironmentModels: getDensity, getDensityBatch!, NoAtmosphereModel
using ..VehicleThermalModels: getHeatRate
using ..ThrusterModels: BaseThrusterModel
using ..DynamicEffectors.AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
using ..GravityEffectors: InverseSquaredJ2GravityModel, j2_secular_rates
using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
using ..ConfigTypes: SaveData
import ..ConfigTypes: GramTrackCache
using ..ControlHooks: calcControlEffect!
using ..GuidanceHooks: calcGuidanceEffect!
using ..NavigationHooks: calcNavigationEffect!
using ..SimConfig: SimulationConfiguration, MissionOrbits
export SaveField, default_save_fields, get_callbacks

const _simulation_model_module = parentmodule(@__MODULE__)

@inline function _simulation_engine_module()
    root = parentmodule(_simulation_model_module)
    isdefined(root, :SimulationEngine) || error("SimulationEngine module is not available for callback stage sampling.")
    return getproperty(root, :SimulationEngine)
end

@inline callback_verbose(integrator) = integrator.p.args.simulation_settings.verbose
@inline callback_use_invokelatest() = get(ENV, "SPACEAGORA_DEV_HOT_RELOAD", "0") == "1"
const _gram_track_cache_warning_emitted = Ref(false)

Base.@kwdef mutable struct GramRuntimeStats
    density_calls::Int64 = 0
    cache_enabled_calls::Int64 = 0
    cache_hits::Int64 = 0
    cache_misses::Int64 = 0
    miss_time_window::Int64 = 0
    miss_state_tolerance::Int64 = 0
    direct_calls::Int64 = 0
    refresh_calls::Int64 = 0
    refresh_points_total::Int64 = 0
    refresh_points_max::Int64 = 0
    refresh_failures::Int64 = 0
    refresh_elapsed_s::Float64 = 0.0
    state_error_samples::Int64 = 0
    alt_err_abs_max_m::Float64 = 0.0
    lat_err_abs_max_deg::Float64 = 0.0
    lon_err_abs_max_deg::Float64 = 0.0
    alt_err_abs_sum_m::Float64 = 0.0
    lat_err_abs_sum_deg::Float64 = 0.0
    lon_err_abs_sum_deg::Float64 = 0.0
end

const _gram_runtime_stats = Ref{GramRuntimeStats}(GramRuntimeStats())
const _gram_runtime_stats_lock = ReentrantLock()

@inline function _gram_runtime_stats_enabled()::Bool
    return _parse_bool_env("SPACEAGORA_GRAM_PROFILE", false)
end

function _gram_runtime_stats_reset!()
    lock(_gram_runtime_stats_lock) do
        _gram_runtime_stats[] = GramRuntimeStats()
    end
    return nothing
end

function _gram_runtime_stats_update!(f::Function)
    lock(_gram_runtime_stats_lock) do
        f(_gram_runtime_stats[])
    end
    return nothing
end

function _gram_runtime_stats_snapshot()
    lock(_gram_runtime_stats_lock) do
        s = _gram_runtime_stats[]
        return (
            density_calls=s.density_calls,
            cache_enabled_calls=s.cache_enabled_calls,
            cache_hits=s.cache_hits,
            cache_misses=s.cache_misses,
            miss_time_window=s.miss_time_window,
            miss_state_tolerance=s.miss_state_tolerance,
            direct_calls=s.direct_calls,
            refresh_calls=s.refresh_calls,
            refresh_points_total=s.refresh_points_total,
            refresh_points_max=s.refresh_points_max,
            refresh_failures=s.refresh_failures,
            refresh_elapsed_s=s.refresh_elapsed_s,
            state_error_samples=s.state_error_samples,
            alt_err_abs_max_m=s.alt_err_abs_max_m,
            lat_err_abs_max_deg=s.lat_err_abs_max_deg,
            lon_err_abs_max_deg=s.lon_err_abs_max_deg,
            alt_err_abs_sum_m=s.alt_err_abs_sum_m,
            lat_err_abs_sum_deg=s.lat_err_abs_sum_deg,
            lon_err_abs_sum_deg=s.lon_err_abs_sum_deg
        )
    end
end
