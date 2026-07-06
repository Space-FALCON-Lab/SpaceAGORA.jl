using Test
using CSV
using DataFrames
using Dates
using LinearAlgebra
using Logging
using Serialization
using StaticArrays
using ComponentArrays
using DiffEqBase
using DiffEqCallbacks
using OrdinaryDiffEq
using Quaternions
using SatelliteToolboxGravityModels
using SPICE
using SpecialFunctions: loggamma
using TOML

const HAS_AQUA = let
    try
        @eval using Aqua
        true
    catch err
        @info "Skipping Aqua-backed test checks" exception=(err, catch_backtrace())
        false
    end
end

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
using .SimulationModel
include(joinpath(REPO_ROOT, "src", "core", "interfaces", "reference_system.jl"))
include(joinpath(REPO_ROOT, "src", "mission", "operations", "maneuver_plans.jl"))

# SimulationEngine uses SimulationModel and provides canonical runtime entrypoints.
const quat_mult = SimulationModel.quat_mult
if !isdefined(@__MODULE__, :SimulationEngine)
    include(joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
end
if !isdefined(@__MODULE__, :run_simulation)
    const run_simulation = SimulationEngine.run_simulation
end
if !isdefined(@__MODULE__, :TelemetryVerification)
    include(joinpath(REPO_ROOT, "src", "analysis", "verification", "telemetry_verification.jl"))
end
if !isdefined(@__MODULE__, :SimulationCampaigns)
    include(joinpath(REPO_ROOT, "src", "simulation", "campaigns", "simulation_campaigns.jl"))
end
if !isdefined(@__MODULE__, :SpaceAGORA)
    include(joinpath(REPO_ROOT, "src", "SpaceAGORA.jl"))
end
const HAS_GRAMSUITE = let
    vendored_gramsuite = joinpath(REPO_ROOT, "data", "GRAMSuite.jl")
    try
        if Base.find_package("GRAMSuite") === nothing && isdir(vendored_gramsuite)
            pushfirst!(LOAD_PATH, vendored_gramsuite)
        end
        @eval import GRAMSuite
        true
    catch err
        @info "Skipping GRAMSuite-backed test checks" exception=(err, catch_backtrace())
        false
    end
end

if HAS_GRAMSUITE
    const TEST_GRAM_LOCK = RuntimeServices.GRAM_LOCK
    const EM = SimulationModel.EnvironmentModels

    EM._GRAM_USE_GLOBAL_LOCK_FN[] = () -> GRAMSuite.gram_use_global_lock()
    EM._GRAM_DEFAULT_SURROGATE_FILE_FN[] = planet -> GRAMSuite.gram_default_surrogate_file(planet)
    EM._CLEAR_GRAM_STATIC_GRID_CACHE_FN[] = () -> GRAMSuite.clear_gram_static_grid_cache!()
    EM._CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[] = () -> GRAMSuite.clear_gram_offline_surrogate_cache!()

    function EM.GRAMAtmosphereModel(; kwargs...)
        return EM.GRAMAtmosphereModel(GRAMSuite.GRAMAtmosphereModel(; kwargs...))
    end

    function EM.GRAMAtmosphereModelSurrogate(;
        surrogate_file::String="",
        point_fallback_below_m::Union{Nothing, Real}=nothing,
        kwargs...
    )
        base_model = EM.GRAMAtmosphereModel(; kwargs...)
        file = isempty(strip(surrogate_file)) ?
            GRAMSuite.gram_default_surrogate_file(base_model.planet_name; gram_root=base_model.gram_root) :
            GRAMSuite.resolve_path(surrogate_file)
        point_fallback = if point_fallback_below_m === nothing
            GRAMSuite.gram_default_point_fallback_below_m(base_model.planet_name)
        else
            value = Float64(point_fallback_below_m)
            value >= 0.0 || throw(ArgumentError("point_fallback_below_m must be >= 0.0 m, got $value"))
            value
        end
        return EM.GRAMAtmosphereModelSurrogate(base_model, file, point_fallback)
    end

    function Base.deepcopy_internal(model::EM.GRAMAtmosphereModel, stackdict::IdDict)
        haskey(stackdict, model) && return stackdict[model]
        copied = lock(TEST_GRAM_LOCK) do
            EM.GRAMAtmosphereModel(deepcopy(model.core))
        end
        stackdict[model] = copied
        return copied
    end

    function Base.deepcopy_internal(model::EM.GRAMAtmosphereModelSurrogate, stackdict::IdDict)
        haskey(stackdict, model) && return stackdict[model]
        copied = lock(TEST_GRAM_LOCK) do
            EM.GRAMAtmosphereModelSurrogate(
                deepcopy(model.base_model),
                model.surrogate_file,
                model.point_fallback_below_m
            )
        end
        stackdict[model] = copied
        return copied
    end

    function EM._gram_core_density_state(
        core::GRAMSuite.GRAMAtmosphereModel,
        h::Float64,
        lat::Float64,
        lon::Float64,
        el_time::Float64,
        wind::Bool,
        lock_obj,
        vacuum_temperature::Float64
    )::Tuple{Float64, Float64, SVector{3, Float64}}
        return GRAMSuite.density_state(
            core, h, lat, lon, el_time, wind;
            lock_obj=lock_obj,
            vacuum_temperature=vacuum_temperature
        )
    end

    @inline function EM._gram_point_density(
        model::EM.GRAMAtmosphereModel,
        h::Float64,
        lat::Float64,
        lon::Float64,
        el_time::Float64,
        wind::Bool
    )::Tuple{Float64, Float64, SVector{3, Float64}}
        h_gram = max(h, -30.0)
        return GRAMSuite.point_density_state(model.core, h_gram, lat, lon, el_time, wind; lock_obj=TEST_GRAM_LOCK)
    end

    function EM.getDensity(
        model::EM.GRAMAtmosphereModel,
        h::Float64,
        lat::Float64,
        lon::Float64,
        el_time::Float64,
        wind::Bool,
        p::params
    )::Tuple{Float64, Float64, SVector{3, Float64}} where {params}
        EI = p.args.environment_model.EI * 1e3
        drag_state = h - EI <= 0.0

        if h > 2000.0e3
            rho = 0.0
            T = p.args.environment_model.planet.T_ref
            wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        elseif !drag_state && !p.args.mission_configuration.keplerian
            rho, T, wind_vec = EM.density_polyfit(h, p)
        else
            h_gram = max(h, -30.0)
            rho, T, wind_vec = GRAMSuite.density_state(
                model.core,
                h_gram,
                lat,
                lon,
                el_time,
                wind;
                lock_obj=TEST_GRAM_LOCK,
                vacuum_temperature=p.args.environment_model.planet.T_ref
            )
        end

        return rho, T, wind_vec
    end

    function EM.getDensity(
        model::EM.GRAMAtmosphereModelSurrogate,
        h::Float64,
        lat::Float64,
        lon::Float64,
        el_time::Float64,
        wind::Bool,
        p::params
    )::Tuple{Float64, Float64, SVector{3, Float64}} where {params}
        EI = p.args.environment_model.EI * 1e3
        drag_state = h - EI <= 0.0

        if h > 2000.0e3
            return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
        elseif !drag_state && !p.args.mission_configuration.keplerian
            return EM.density_polyfit(h, p)
        end

        base_model = model.base_model isa EM.GRAMAtmosphereModel ? model.base_model.core : model.base_model
        point_fallback = model.base_model isa EM.GRAMAtmosphereModel ? nothing :
            (m, h_i, lat_i, lon_i, t_i, w_i) -> EM._gram_point_density(m, h_i, lat_i, lon_i, t_i, w_i)
        h_gram = max(h, -30.0)

        return GRAMSuite.surrogate_density_state(
            base_model,
            model.surrogate_file,
            model.point_fallback_below_m,
            h_gram,
            lat,
            lon,
            el_time,
            wind;
            lock_obj=TEST_GRAM_LOCK,
            point_density_fallback=point_fallback,
            vacuum_temperature=p.args.environment_model.planet.T_ref
        )
    end
end

if !isdefined(@__MODULE__, :make_example_config)
    const make_example_config = TelemetryVerification.make_example_config
    const make_three_body_spacecraft = TelemetryVerification.make_three_body_spacecraft
    const run_and_report = TelemetryVerification.run_and_report
end
if !isdefined(@__MODULE__, :_solver_policy_mode)
    const build_initial_conditions = SimulationEngine.build_initial_conditions
    const _build_solver_tolerances = SimulationEngine._build_solver_tolerances
    const _solve_with_solver_policy = SimulationEngine._solve_with_solver_policy
    const _split_subproblem = SimulationEngine._split_subproblem
    const _resolve_component_tolerance = SimulationEngine._resolve_component_tolerance
    const spacecraft_dynamics! = SimulationEngine.spacecraft_dynamics!
    const _debug_print_nan_parameter_paths! = SimulationEngine._debug_print_nan_parameter_paths!
    const _load_checkpoint = SimulationEngine._load_checkpoint
    const _normalize_warning_emitted = SimulationEngine._normalize_warning_emitted
    const _solver_policy_mode = SimulationEngine._solver_policy_mode
    const _gram_per_sat_instances_enabled = SimulationEngine._gram_per_sat_instances_enabled
    const _solver_maxiters = SimulationEngine._solver_maxiters
    const _effector_parallel_mode = SimulationEngine._effector_parallel_mode
    const _effector_thread_threshold = SimulationEngine._effector_thread_threshold
    const _effector_max_threads = SimulationEngine._effector_max_threads
    const _effector_long_mission_threshold_s = SimulationEngine._effector_long_mission_threshold_s
    const _effector_cost_ns_per_item_default = SimulationEngine._effector_cost_ns_per_item_default
    const _effector_cost_ema_alpha = SimulationEngine._effector_cost_ema_alpha
    const _effector_work_ns_per_worker_threshold = SimulationEngine._effector_work_ns_per_worker_threshold
    const _dynamic_effector_thread_decision = SimulationEngine._dynamic_effector_thread_decision
    const _rhs_execution_mode_env = SimulationEngine._rhs_execution_mode_env
    const _rhs_execution_plan = SimulationEngine._rhs_execution_plan
    const _retcode_is_stiff_symptom = SimulationEngine._retcode_is_stiff_symptom
    const _split_imex_solver_spec = SimulationEngine._split_imex_solver_spec
    const _symplectic_conservative_eligible = SimulationEngine._symplectic_conservative_eligible
    const _symplectic_fixed_dt_s = SimulationEngine._symplectic_fixed_dt_s
    const _gravity_backbone_fixed_dt_s = SimulationEngine._gravity_backbone_fixed_dt_s
    const _gravity_backbone_eligible = SimulationEngine._gravity_backbone_eligible
    const _gravity_backbone_reject_reason = SimulationEngine._gravity_backbone_reject_reason
    const _gravity_backbone_structure_validated = SimulationEngine._gravity_backbone_structure_validated
    const _gravity_backbone_kick_structure_validated = SimulationEngine._gravity_backbone_kick_structure_validated
    const _auto_stiff_switched = SimulationEngine._auto_stiff_switched
    const _solve_with_explicit_solver = SimulationEngine._solve_with_explicit_solver
    const _solve_with_multirate_solver = SimulationEngine._solve_with_multirate_solver
    const _multirate_fast_substeps = SimulationEngine._multirate_fast_substeps
    const _multirate_slow_dt_s = SimulationEngine._multirate_slow_dt_s
    const _multirate_slow_solver_spec = SimulationEngine._multirate_slow_solver_spec
    const _multirate_fast_solver_spec = SimulationEngine._multirate_fast_solver_spec
    const _append_series_columns! = SimulationEngine._append_series_columns!
    const _atomic_write_file = SimulationEngine._atomic_write_file
    const _checkpoint_directory = SimulationEngine._checkpoint_directory
    const _checkpoint_paths = SimulationEngine._checkpoint_paths
    const _clear_ephemeris_reuse_cache! = SimulationEngine._clear_ephemeris_reuse_cache!
    const _effector_long_orbit_threshold = SimulationEngine._effector_long_orbit_threshold
    const _ephemeris_reuse_enabled = SimulationEngine._ephemeris_reuse_enabled
    const _ephemeris_reuse_max_entries = SimulationEngine._ephemeris_reuse_max_entries
    const _ephemeris_reuse_store! = SimulationEngine._ephemeris_reuse_store!
    const _find_sample_value = SimulationEngine._find_sample_value
    const _has_active_srp_effector = SimulationEngine._has_active_srp_effector
    const _initialize_aero_workspace_buffers! = SimulationEngine._initialize_aero_workspace_buffers!
    const _initialize_harmonics_workspace_buffers! = SimulationEngine._initialize_harmonics_workspace_buffers!
    const _initialize_heat_rate_buffers! = SimulationEngine._initialize_heat_rate_buffers!
    const _initialize_nbody_ephemeris_cache! = SimulationEngine._initialize_nbody_ephemeris_cache!
    const _initialize_nbody_ephemeris_cache_buffer! = SimulationEngine._initialize_nbody_ephemeris_cache_buffer!
    const _initialize_nbody_workspace_buffers! = SimulationEngine._initialize_nbody_workspace_buffers!
    const _initialize_planet_frame_cache_buffer! = SimulationEngine._initialize_planet_frame_cache_buffer!
    const _initialize_planet_frame_ephemeris_cache! = SimulationEngine._initialize_planet_frame_ephemeris_cache!
    const _initialize_srp_sun_cache_buffer! = SimulationEngine._initialize_srp_sun_cache_buffer!
    const _initialize_srp_sun_ephemeris_cache! = SimulationEngine._initialize_srp_sun_ephemeris_cache!
    const _initialize_spice_rhs_memo_mode! = SimulationEngine._initialize_spice_rhs_memo_mode!
    const _mission_is_long_for_effector_threads = SimulationEngine._mission_is_long_for_effector_threads
    const _nbody_ephemeris_cache_dt_s = SimulationEngine._nbody_ephemeris_cache_dt_s
    const _nbody_ephemeris_cache_max_samples = SimulationEngine._nbody_ephemeris_cache_max_samples
    const _planet_frame_cache_dt_s = SimulationEngine._planet_frame_cache_dt_s
    const _planet_frame_cache_max_samples = SimulationEngine._planet_frame_cache_max_samples
    const _reset_spice_rhs_memo! = SimulationEngine._reset_spice_rhs_memo!
    const _reset_spice_runtime_counters! = SimulationEngine._reset_spice_runtime_counters!
    const _srp_ephemeris_cache_dt_s = SimulationEngine._srp_ephemeris_cache_dt_s
    const _srp_ephemeris_cache_max_samples = SimulationEngine._srp_ephemeris_cache_max_samples
    const _validate_orientation_inertia! = SimulationEngine._validate_orientation_inertia!
    const _validate_thermal_model_support! = SimulationEngine._validate_thermal_model_support!
    const build_state_sample = SimulationEngine.build_state_sample
    const sample_planet_frame = SimulationEngine.sample_planet_frame
    const sample_atmosphere = SimulationEngine.sample_atmosphere
    const sample_environment = SimulationEngine.sample_environment
    const _build_typed_solver_problem = SimulationEngine._build_typed_solver_problem
    const _evaluate_dynamic_effector = SimulationEngine._evaluate_dynamic_effector
    const _solver_partition_validated = SimulationEngine._solver_partition_validated
    const _wrench_method_available = SimulationEngine._wrench_method_available
    const spacecraft_dynamics_gravity_backbone! = SimulationEngine.spacecraft_dynamics_gravity_backbone!
    const _gravity_backbone_half_kick! = SimulationEngine._gravity_backbone_half_kick!
    const spacecraft_dynamics_implicit_atmosphere! = SimulationEngine.spacecraft_dynamics_implicit_atmosphere!
    const spacecraft_dynamics_explicit_remainder! = SimulationEngine.spacecraft_dynamics_explicit_remainder!
    const _is_gravity_backbone_state = SimulationEngine._is_gravity_backbone_state
    const _state_position_ii = SimulationEngine._state_position_ii
    const _state_velocity_ii = SimulationEngine._state_velocity_ii
    const _state_mass_kg = SimulationEngine._state_mass_kg
    const _state_heat_loads = SimulationEngine._state_heat_loads
    const _gravity_backbone_initial_states = SimulationEngine._gravity_backbone_initial_states
end

struct ConstantDensityModel <: SimulationModel.AbstractDensityModel
    rho::Float64
    temp::Float64
end

struct TimedTangentialThrusterModel <: SimulationModel.AbstractControlEffectorModel
    thrust::Float64
    direction_sign::Float64 # +1.0 prograde, -1.0 retrograde
    start_time::Float64
    stop_time::Float64
end

mutable struct CountingGuidanceModel <: SimulationModel.AbstractTypes.AbstractGuidanceModel
    hits::Vector{Int}
end

mutable struct CountingNavigationModel
    hits::Vector{Int}
end

mutable struct CountingControlModel <: SimulationModel.AbstractControlEffectorModel
    hits::Vector{Int}
end

struct ThrowingForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct NaNForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct NaNParamForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct ConstantTorqueModel <: SimulationModel.AbstractForceTorqueModel
    torque::SVector{3, Float64}
end

struct ConstantForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
end

struct WrenchOnlyForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
    torque::SVector{3, Float64}
end

struct AtmosphereProbeWrenchModel <: SimulationModel.AbstractForceTorqueModel
end

struct ImplicitLegacyForceModel <: SimulationModel.AbstractForceTorqueModel
    force::SVector{3, Float64}
    torque::SVector{3, Float64}
end

struct InvalidPartitionForceModel <: SimulationModel.AbstractForceTorqueModel
end

struct BackboneCustomGravityModel <: SimulationModel.AbstractForceTorqueModel
    accel::SVector{3, Float64}
end

struct InvalidBackboneStructureModel <: SimulationModel.AbstractForceTorqueModel
end

struct SolarBackboneModel <: SimulationModel.AbstractForceTorqueModel
end

struct InvalidBackboneKickStructureModel <: SimulationModel.AbstractForceTorqueModel
end

struct PlanetFrameKickModel <: SimulationModel.AbstractForceTorqueModel
end

struct ThrowingOrbitPlanet <: SimulationModel.AbstractPlanet
    Rp_e::Float64
end

function SimulationModel.EnvironmentModels.getDensity(
    model::ConstantDensityModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)
    return model.rho, model.temp, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlForceTorque(
    model::TimedTangentialThrusterModel,
    u::AbstractVector{Float64},
    p::ODEParams,
    i::Int64,
    t::Float64
)
    if t < model.start_time || t > model.stop_time
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    v = SVector{3, Float64}(u.vel)
    vm = norm(v)
    if vm == 0.0 || !isfinite(vm)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    force = model.thrust * model.direction_sign * (v / vm)
    return force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcControlEffect!(
    model::TimedTangentialThrusterModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    return nothing
end

function SimulationModel.calcGuidanceEffect!(
    model::CountingGuidanceModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcNavigationEffect!(
    model::CountingNavigationModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcControlEffect!(
    model::CountingControlModel,
    u::ComponentVector,
    p::ODEParams,
    t::Float64,
    i::Int64
)
    model.hits[i] += 1
    return nothing
end

function SimulationModel.calcForceTorque(
    model::ThrowingForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    error("intentional derivative failure")
end

function SimulationModel.calcForceTorque(
    model::NaNForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(NaN, NaN, NaN), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::NaNParamForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    p.shared_buffers.current_time[] = NaN
    return SVector{3, Float64}(NaN, NaN, NaN), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::ConstantTorqueModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), model.torque
end

function SimulationModel.calcForceTorque(
    model::ConstantForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return model.force, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.wrench(
    model::WrenchOnlyForceModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return model.force, model.torque
end

SimulationModel.solver_partition(::ImplicitLegacyForceModel) = :implicit
SimulationModel.solver_partition(::InvalidPartitionForceModel) = :bad_partition
SimulationModel.gravity_backbone_structure(::BackboneCustomGravityModel) = :position_only_static_gravity
SimulationModel.gravity_backbone_structure(::InvalidBackboneStructureModel) = :bad_structure
SimulationModel.gravity_backbone_structure(::SolarBackboneModel) = :position_only_static_gravity
SimulationModel.gravity_backbone_kick_structure(::InvalidBackboneKickStructureModel) = :bad_structure
SimulationModel.gravity_backbone_kick_structure(::PlanetFrameKickModel) = :velocity_kick_explicit

function SimulationModel.gravity_backbone_acceleration_ii(
    model::BackboneCustomGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return model.accel
end

function SimulationModel.gravity_backbone_acceleration_ii(
    model::SolarBackboneModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.environment_requirements(::SolarBackboneModel)
    return EffectorEnvironmentRequirements(solar=true)
end

function SimulationModel.environment_requirements(::PlanetFrameKickModel)
    return EffectorEnvironmentRequirements(planet_frame=true)
end

function SimulationModel.gravity_backbone_kick_acceleration_ii(
    model::PlanetFrameKickModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.calcForceTorque(
    model::ImplicitLegacyForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return model.force, model.torque
end

function SimulationModel.calcForceTorque(
    model::InvalidPartitionForceModel,
    x::AbstractVector{Float64},
    p::ODEParams,
    i::Int64
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.environment_requirements(::AtmosphereProbeWrenchModel)
    return EffectorEnvironmentRequirements(planet_frame=true, atmosphere=true)
end

function SimulationModel.wrench(
    model::AtmosphereProbeWrenchModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64
)
    env.atmosphere === nothing && error("expected atmosphere sample")
    env.planet_frame === nothing && error("expected planet-frame sample")
    return SVector{3, Float64}(env.atmosphere.rho_kg_m3, env.planet_frame.alt_m, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function SimulationModel.ControlHooks.rvtoorbitalelement(
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    planet::ThrowingOrbitPlanet
)
    throw(ErrorException("forced-orbital-element-conversion-failure"))
end

const SPICE_PATH = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
const EARTH = Earth("", SPICE_PATH)

angle_distance(a::Float64, b::Float64) = abs(mod(a - b + pi, 2pi) - pi)

function specific_energy(df::DataFrame, mu::Float64)
    r = sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    v2 = df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2
    return 0.5 .* v2 .- mu ./ r
end

function specific_angular_momentum_magnitude(df::DataFrame)
    h1 = df.sc1_pos_2 .* df.sc1_vel_3 .- df.sc1_pos_3 .* df.sc1_vel_2
    h2 = df.sc1_pos_3 .* df.sc1_vel_1 .- df.sc1_pos_1 .* df.sc1_vel_3
    h3 = df.sc1_pos_1 .* df.sc1_vel_2 .- df.sc1_pos_2 .* df.sc1_vel_1
    return sqrt.(h1.^2 .+ h2.^2 .+ h3.^2)
end

function unwrap_angle_series(values::AbstractVector{<:Real})
    n = length(values)
    n == 0 && return Float64[]
    out = Vector{Float64}(undef, n)
    out[1] = Float64(values[1])
    @inbounds for i in 2:n
        prev = Float64(values[i - 1])
        cur = Float64(values[i])
        Δ = atan(sin(cur - prev), cos(cur - prev))
        out[i] = out[i - 1] + Δ
    end
    return out
end

function linear_regression_slope(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    length(xs) == length(ys) || throw(ArgumentError("linear_regression_slope requires vectors of equal length."))
    length(xs) >= 2 || throw(ArgumentError("linear_regression_slope requires at least two samples."))
    x̄ = sum(xs) / length(xs)
    ȳ = sum(ys) / length(ys)
    num = 0.0
    den = 0.0
    @inbounds for i in eachindex(xs, ys)
        dx = Float64(xs[i]) - x̄
        dy = Float64(ys[i]) - ȳ
        num += dx * dy
        den += dx * dx
    end
    den > 0.0 || throw(ArgumentError("linear_regression_slope requires non-degenerate x samples."))
    return num / den
end

function harmonics_full_normalization_factor(l::Int, m::Int)::Float64
    l >= 0 || throw(ArgumentError("Degree must be >= 0, got $l."))
    0 <= m <= l || throw(ArgumentError("Order must satisfy 0 <= m <= l, got l=$l, m=$m."))
    δ0m = (m == 0) ? 1.0 : 0.0
    ratio = exp(0.5 * (loggamma(l - m + 1) - loggamma(l + m + 1)))
    return sqrt((2.0 - δ0m) * (2.0 * l + 1.0)) * ratio
end

function write_gravity_harmonics_csv(
    coeff_fn::Function,
    path::AbstractString,
    max_degree::Int,
    max_order::Int
)
    open(path, "w") do io
        println(io, "degree,order,C,S")
        for l in 0:max_degree
            for m in 0:min(max_order, l)
                C, S = coeff_fn(l, m)
                println(io, string(l, ",", m, ",", C, ",", S))
            end
        end
    end
    return path
end

function write_icgem_from_harmonics_model(
    path::AbstractString,
    model::SimulationModel.DynamicEffectors.GravitationalHarmonicsModel;
    model_name::AbstractString="SpaceAGORA_Test"
)
    open(path, "w") do io
        println(io, "product_type gravity_field")
        println(io, "modelname ", model_name)
        println(io, "earth_gravity_constant ", model.planet.μ)
        println(io, "radius ", model.planet.Rp_e)
        println(io, "max_degree ", model.L)
        println(io, "errors no")
        println(io, "norm fully_normalized")
        println(io, "end_of_head")
        for l in 0:model.L
            for m in 0:min(model.M, l)
                i = l + 1
                j = m + 1
                println(io, "gfc ", l, " ", m, " ", model.C[i, j], " ", model.S[i, j])
            end
        end
    end
    return path
end

function make_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0,
    orientation_state::Union{Nothing, Tuple{SVector{4, Float64}, SVector{3, Float64}}}=nothing
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    panel = Link{0}(root=false, m=30.0, ref_area=6.0, r=MVector{3, Float64}(0.0, 1.2, 0.0))

    if isnothing(orientation_state)
        ic = InitialCondition(
            ra=EARTH.Rp_e + ra_alt_m,
            rp=EARTH.Rp_e + rp_alt_m,
            i=i_deg,
            ω=ω_deg,
            Ω=Ω_deg,
            ν=ν_deg
        )
    else
        q0, w0 = orientation_state
        ra = EARTH.Rp_e + ra_alt_m
        rp = EARTH.Rp_e + rp_alt_m
        a = (ra + rp) / 2.0
        e = (ra - rp) / (ra + rp)
        ic = InitialCondition(a, e, i_deg, ω_deg, Ω_deg, ν_deg, q0, w0)
    end

    dry_mass = root.m + panel.m
    return SpacecraftModel(
        Joint[],
        [root, panel],
        root,
        true,
        dry_mass,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function make_single_link_spacecraft(;
    ra_alt_m::Float64,
    rp_alt_m::Float64,
    i_deg::Float64=35.0,
    ω_deg::Float64=40.0,
    Ω_deg::Float64=10.0,
    ν_deg::Float64=175.0
)
    root = Link{0}(root=true, m=500.0, ref_area=12.0)
    ic = InitialCondition(
        ra=EARTH.Rp_e + ra_alt_m,
        rp=EARTH.Rp_e + rp_alt_m,
        i=i_deg,
        ω=ω_deg,
        Ω=Ω_deg,
        ν=ν_deg
    )

    return SpacecraftModel(
        Joint[],
        [root],
        root,
        true,
        root.m,
        0.0,
        root.inertia,
        0,
        0,
        ic,
        1
    )
end

function make_base_thruster_model(;
    thrust::Real,
    direction::Real=0.0,
    Δv::Real=0.0,
    start_burn_time::Real=0.0,
    stop_burn_time::Real=0.0,
    Isp::Real=300.0
)
    return BaseThrusterModel(
        thrust=[Float64(thrust)],
        direction=[Float64(direction)],
        Δv=[Float64(Δv)],
        start_burn_time=[Float64(start_burn_time)],
        stop_burn_time=[Float64(stop_burn_time)],
        Isp=[Float64(Isp)]
    )
end

function build_config(;
    spacecraft::SpacecraftModel,
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::SimulationModel.InitialTime=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    ephemerides_model=SpiceEphemeridesModel(),
    planet=EARTH
)
    environment_model = EnvironmentModel(
        planet=planet,
        EI=EI_km,
        density_model=density_model,
        ephemerides_model=ephemerides_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    )

    return SimulationConfiguration(
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel([spacecraft], dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

function build_config_multi(;
    spacecraft::Vector{SpacecraftModel},
    density_model,
    orientation_sim::Bool,
    mission_time::Float64,
    EI_km::Float64,
    dynamic_effectors::Tuple,
    control_effectors::Tuple=(),
    control_rates::Vector{Float64}=Float64[],
    keplerian::Bool=true,
    simulation_settings::SimulationSettings=SimulationSettings(results=true, verbose=false, generate_plots=false, normalize=false),
    tolerances::IntegrationTolerances=IntegrationTolerances(),
    initial_time::SimulationModel.InitialTime=SimulationModel.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
    ephemerides_model=SpiceEphemeridesModel(),
    planet=EARTH
)
    environment_model = EnvironmentModel(
        planet=planet,
        EI=EI_km,
        density_model=density_model,
        ephemerides_model=ephemerides_model,
        thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
        topography=false,
        wind=false
    )

    return SimulationConfiguration(
        simulation_settings=simulation_settings,
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=keplerian,
            number_of_orbits=1,
            mission_time=mission_time,
            orientation_sim=orientation_sim,
            num_steps_to_save=1000
        ),
        environment_model=environment_model,
        dynamics_model=DynamicsModel(spacecraft, dynamic_effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=control_effectors, control_rates=control_rates),
        initial_time=initial_time,
        integration_tolerances=tolerances
    )
end

function interp_linear(times::AbstractVector{<:Real}, values::AbstractVector{<:Real}, t::Float64)
    idx = searchsortedlast(times, t)
    if idx <= 0
        return Float64(values[1])
    elseif idx >= length(times)
        return Float64(values[end])
    end
    t1 = Float64(times[idx])
    t2 = Float64(times[idx + 1])
    y1 = Float64(values[idx])
    y2 = Float64(values[idx + 1])
    α = (t - t1) / (t2 - t1)
    return (1.0 - α) * y1 + α * y2
end

function make_agora_earth_spacecraft()
    main_bus = Link{0}(root=true, m=620.0, ref_area=2.05 * 2.8)
    left_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, -2.05 / 2.0 - 5.7 / 4.0, 0.0))
    right_panel = Link{0}(root=false, m=10.0, ref_area=5.7 * 1.0 / 2.0, r=MVector{3, Float64}(0.0, 2.05 / 2.0 + 5.7 / 4.0, 0.0))
    ic = InitialCondition(
        ra=56_378.7978559e3,
        rp=EARTH.Rp_e + 200_590.0,
        i=89.876,
        ω=75.505,
        Ω=104.115,
        ν=175.0
    )

    return SpacecraftModel(
        Joint[],
        [main_bus, left_panel, right_panel],
        main_bus,
        true,
        main_bus.m + left_panel.m + right_panel.m,
        200.0,
        main_bus.inertia,
        0,
        0,
        ic,
        1
    )
end

function run_case(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            run_simulation(args; isolate_state=isolate_state)
            @test isfile(results_csv_path)
            return CSV.read(results_csv_path, DataFrame)
        end
    end
end

function run_case_silent(args::SimulationConfiguration; isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            redirect_stdout(devnull) do
                run_simulation(args; isolate_state=isolate_state)
            end
            @test isfile(results_csv_path)
            return CSV.read(results_csv_path, DataFrame)
        end
    end
end

function run_case_capture_stdout(args::SimulationConfiguration; expect_results_csv::Bool=true, isolate_state::Bool=true)
    return mktempdir() do tmp
        cd(tmp) do
            results_csv_path = joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
            output = ""
            mktemp() do path, io
                redirect_stdout(io) do
                    run_simulation(args; isolate_state=isolate_state)
                end
                flush(io)
                seekstart(io)
                output = read(io, String)
            end
            if expect_results_csv
                @test isfile(results_csv_path)
                return CSV.read(results_csv_path, DataFrame), output
            else
                @test !isfile(results_csv_path)
                return DataFrame(), output
            end
        end
    end
end

function seed_solution_for_save_csv!(
    solution::Solution;
    n_bodies::Int,
    n_reaction_wheels::Int,
    n_thrusters::Int,
    base::Float64=1.0,
    closed_form::Bool=false
)
    solution.physical_properties.α = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.β = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_rate = [Float64[] for _ in 1:n_bodies]
    solution.performance.heat_load = [Float64[] for _ in 1:n_bodies]
    solution.physical_properties.rw_h = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.rw_τ = [Float64[] for _ in 1:n_reaction_wheels]
    solution.physical_properties.thruster_forces = [Float64[] for _ in 1:n_thrusters]

    function _push_sample!(field)
        if field isa Vector{Float64}
            push!(field, base)
        elseif field isa Vector{Int64}
            push!(field, 1)
        elseif field isa Vector{Vector{Float64}}
            for subfield in field
                push!(subfield, base)
            end
        end
        return nothing
    end

    for group in (solution.orientation, solution.physical_properties, solution.performance, solution.forces)
        for fname in fieldnames(typeof(group))
            _push_sample!(getfield(group, fname))
        end
    end

    if closed_form
        solution.closed_form.t_cf = [base + 10.0]
        solution.closed_form.h_cf = [base + 20.0]
        solution.closed_form.γ_cf = [base + 30.0]
        solution.closed_form.v_cf = [base + 40.0]
    end

    return solution
end


module IncludeOrderSandbox
end
const INCLUDE_ORDER_SANDBOX = IncludeOrderSandbox

module ExportImportSandbox
end
const EXPORT_IMPORT_SANDBOX = ExportImportSandbox

const GUIDANCE_SANDBOX_LOADED = Ref(false)
module GuidanceSandbox
using ..SimulationModel
using ..SimulationModel.AbstractTypes: AbstractGuidanceModel
using ComponentArrays
using LinearAlgebra
using StaticArrays
end
const GUIDANCE_SANDBOX = GuidanceSandbox


include(joinpath(REPO_ROOT, "test", "suites", "01_contract_and_api_tests.jl"))
include(joinpath(REPO_ROOT, "test", "suites", "02_callbacks_parallel_and_smoke_tests.jl"))
include(joinpath(REPO_ROOT, "test", "suites", "03_persistence_units_and_rotational_tests.jl"))
include(joinpath(REPO_ROOT, "test", "suites", "04_solver_env_and_regression_tests.jl"))
include(joinpath(REPO_ROOT, "test", "suites", "05_thruster_control_and_quality_tests.jl"))
include(joinpath(REPO_ROOT, "test", "suites", "06_monolith_split_runtime_tests.jl"))
include(joinpath(REPO_ROOT, "test", "suites", "07_no_gram_onboarding_tests.jl"))
include(joinpath(REPO_ROOT, "test", "suites", "08_cli_and_assets_tests.jl"))
