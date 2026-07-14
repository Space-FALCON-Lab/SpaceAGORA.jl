# Guards against double-inclusion: golden and legacy-suite test entrypoints
# both need this bootstrap, and re-including it in the same session would
# redefine structs/modules (ConstantDensityModel, GuidanceSandbox, ...) that
# have no per-file isdefined guard of their own.
if !isdefined(Main, :SPACEAGORA_TEST_BOOTSTRAP_LOADED)
const SPACEAGORA_TEST_BOOTSTRAP_LOADED = true

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
if !isdefined(@__MODULE__, :ParallelProfiles)
    # SimulationCampaigns consumes the outer-route bandit via `..ParallelProfiles`.
    include(joinpath(REPO_ROOT, "src", "parallel", "routing", "parallel_profiles.jl"))
end
if !isdefined(@__MODULE__, :ParallelProcess)
    # SimulationCampaigns consumes the process worker pool via `..ParallelProcess`.
    include(joinpath(REPO_ROOT, "src", "parallel", "process", "parallel_process.jl"))
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
    const _active_solver_config = SimulationEngine._active_solver_config
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

include(joinpath(REPO_ROOT, "test", "helpers", "mock_models.jl"))
include(joinpath(REPO_ROOT, "test", "helpers", "test_constants.jl"))
include(joinpath(REPO_ROOT, "test", "helpers", "numeric_helpers.jl"))
include(joinpath(REPO_ROOT, "test", "helpers", "spacecraft_builders.jl"))
include(joinpath(REPO_ROOT, "test", "helpers", "config_builders.jl"))
include(joinpath(REPO_ROOT, "test", "helpers", "simulation_run_helpers.jl"))
include(joinpath(REPO_ROOT, "test", "helpers", "persistence_fixtures.jl"))
include(joinpath(REPO_ROOT, "test", "helpers", "sandbox_modules.jl"))

end # SPACEAGORA_TEST_BOOTSTRAP_LOADED
