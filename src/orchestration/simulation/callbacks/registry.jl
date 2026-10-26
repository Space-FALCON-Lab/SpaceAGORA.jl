using ..FrameTransforms
import ...SimulationLifecycle
import ...GNCRuntime
using ...GNCRuntime: tick_schedule, gnc_callbacks, register_tstops!, gnc_clock_started!, gnc_period_s
using ...GNCRuntime: Produces, Consumes, gnc_links, link_cardinality, check_link, link_initialize!

using OrdinaryDiffEq
using DiffEqBase
using DiffEqCallbacks
using LinearAlgebra
using StaticArrays
using SPICE
using Dates
using ...RuntimeServices: SPICE_LOCK, GRAM_LOCK, tracked_lock
using ..SimulationModel: PlanetFrameEphemerisCache, rot
using ..SimulationModel: ExternalPropagation
using ..ArticulatedBody: articulated_forward_kinematics, articulated_link_poses, articulated_has_moving_joints
using ..SimulationModel: ephemerides_requires_spice
using ..InterLinkModels: has_ground_stations, has_station_counters
using ..SimulationModel.EphemeridesModels: ephemerides_epoch_seconds, ephemerides_frame_lpi, ephemerides_sun_ii
using ..ParallelPolicy
using ..EnvironmentModels
using ..EnvironmentModels: getDensity, getDensityBatch!, NoAtmosphereModel
using ..VehicleThermalModels: getHeatRate
using ..ThrusterModels: BaseThrusterModel
using ..DynamicEffectors: PlumeSurfaceInteraction
using ..DynamicEffectors.AerodynamicEffectors: AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
using ..GravityEffectors: InverseSquaredJ2GravityModel, j2_secular_rates
using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
using ..ConfigTypes: SaveData
using ..SceneVisualization: link_pose_vector, link_pose_link_indices, arm_pose_vector, robot_arm_plan_for
import ..ConfigTypes: GramTrackCache, VacuumPredictedGRAMCache, GramDensityPerturbationState
import ..ConfigTypes: GramTrackCacheConfig, CallbackEnvConfig, PolicyDecisionEnvConfig
using ..ControlHooks: update_control!, control_thruster_levels, _thruster_levels_at, touchdown_spec, refresh_burn!, mounted_burn_tick!, burn_edge
using ..ControlHooks: _effective_burn_window, _burn_window_changed, _burn_window_pending
using ..Actuators: actuation_mode, actuator_levels
using ..ActuatorMounting: mounted_gimbal_count, apply_gimbal_levels!
using ..Components: _gimbal_rotation
using ..ActuatorMounting: mounted_thruster_fractions, mounted_hinge_count, apply_hinge_levels!, mounted_actuator_count
using ..ActuatorMounting: applied_actuator_levels, mounted_power
using ..ActuatorMounting: _has_run_state_parts, _note_run_writes!, _start_run_state!
using ..Actuators: actuator_saturate
using ..ThrusterHooks: thrust_calculation_schmitt_trigger!
using ..PowerSystems: has_battery, has_power_parts, battery, solar_generation_w, solar_flux_w_m2, mean_solar_distance_m
using ..PowerSystems: bus_load_w, actuator_power_by_tier, battery_boundary!, update_power_state!, shed_levels, within_full_tolerance
using ..InterlinkTerminals: without_link_driven, link_driven_levels!, link_terminal_power, link_driven
using ..InterlinkTerminals: link_receives, received_link_power_w, beamed_link_power_w
using ..InterlinkTerminals: data_link_rates_bps, shed_link_driven_levels!, _data_flows
using ..DataStorage: has_data_store, data_store, _data_store_boundary!, _update_data_state!, _start_data_state!
using ..InterLinkModels: LinkContext, PowerBeamLink, has_data_links
using ..DynamicEffectors: eclipse_area_calc, MagneticTorqueRodModel, rod_field_body
using ..DynamicEffectors: is_gravitational
using ..ActuatorFaults: has_faults, faulted_levels, refresh_faults!, fault_edge, resolve_faults!
using ..Sensors: SensorTruth, sensor_measure, sensor_requirements, sensor_period_s, initialize_sensor!, sensor_needs_specific_force
using ..SensorMounting: SensorMounting, mounted_sensors, reset_sensors!
using ..ActuatorMounting: actuator_context
using ..TerrainModels: AbstractTerrainModel, DEMTerrainModel, terrain_radius
using ..ControlHooks: _edg_heat_states, _edg_heat_boundary_distance, _edg_heat_boundary_vertical_speed, _edg_capture_entry_heat!, _edg_capture_exit_heat!, _edg_reconcile_heat_boundary!
using ..GuidanceHooks: update_guidance!
using ..NavigationHooks: update_navigation!
using ..SimConfig: SimulationConfiguration, MissionOrbits
export SaveField, default_save_fields, get_callbacks
export TrajectoryRecorder, get_trajectory_recorder_callback
export record_trajectory_sample!, reset_trajectory_recorder!
export trajectory_times, trajectory_field, trajectory_save_data
export trajectory_positions, trajectory_velocities, trajectory_masses

const _simulation_model_module = parentmodule(@__MODULE__)
import ...CallbackServices

@inline callback_verbose(integrator) = integrator.p.args.simulation_settings.verbose
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
