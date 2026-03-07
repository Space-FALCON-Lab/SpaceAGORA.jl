include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
using LinearAlgebra
using StaticArrays
using LoopVectorization
using ComponentArrays
using DifferentialEquations
using CSV
using DataFrames
using Polyester
using Arrow
using Dates
using Serialization
using SHA
using TOML

const _normalize_warning_emitted = Ref(false)
const RESULTS_BUNDLE_SCHEMA_VERSION = "1"
const CHECKPOINT_SCHEMA_VERSION = "1"
const _EPHEMERIS_REUSE_LOCK = ReentrantLock()
const _SRP_EPHEMERIS_REUSE_CACHE = Dict{Any, SimulationModel.SRPSunEphemerisCache}()
const _NBODY_EPHEMERIS_REUSE_CACHE = Dict{Any, SimulationModel.NBodyEphemerisCache}()
const _PLANET_FRAME_EPHEMERIS_REUSE_CACHE = Dict{Any, SimulationModel.PlanetFrameEphemerisCache}()

@inline _typed_normalize_warning_enabled() = _engine_env_get("SPACEAGORA_WARN_NORMALIZE", "1") == "1"
@inline _typed_allow_legacy_normalize() = _engine_env_get("SPACEAGORA_ALLOW_TYPED_NORMALIZE", "0") == "1"
@inline _typed_save_bundle_enabled() = _engine_env_get("SPACEAGORA_SAVE_BUNDLE", "1") == "1"
@inline _typed_checkpoint_enabled(args) = args.simulation_settings.checkpoint_enabled || args.simulation_settings.resume_from_checkpoint

function _validate_orientation_inertia!(args)
    if !args.mission_configuration.orientation_sim
        return nothing
    end
    for (i, sc) in enumerate(args.dynamics_model.spacecraft)
        inertia_tensor = sc.inertia_tensor
        inertia_matrix = Matrix(inertia_tensor)
        if !all(isfinite, inertia_matrix)
            throw(ArgumentError("Spacecraft index $i has non-finite inertia tensor entries, required for orientation_sim=true."))
        end
        if !issymmetric(inertia_matrix)
            throw(ArgumentError("Spacecraft index $i has non-symmetric inertia tensor, required for orientation_sim=true."))
        end
        if !isposdef(Symmetric(inertia_matrix))
            throw(ArgumentError("Spacecraft index $i has non-positive-definite inertia tensor, required for orientation_sim=true."))
        end
    end
    return nothing
end

function _validate_thermal_model_support!(args)
    thermal_model = args.environment_model.thermal_model
    if !hasmethod(SimulationModel.getHeatRate, Tuple{typeof(thermal_model), Float64, Float64, Float64, Float64, Float64})
        throw(ArgumentError(
            "Thermal model $(nameof(typeof(thermal_model))) must implement " *
            "getHeatRate(model, S, T, ρ, v, α)::Float64."
        ))
    end
    return nothing
end

function _validate_ephemerides_support!(args)
    ephemerides_model = args.environment_model.ephemerides_model
    if ephemerides_model isa SimulationModel.SimpleEphemeridesModel
        if _has_active_nbody_effector(args.dynamics_model.dynamic_effectors)
            throw(ArgumentError(
                "SimpleEphemeridesModel does not support active N-body gravity effectors. " *
                "Use SpiceEphemeridesModel() for high-fidelity third-body runs."
            ))
        end
        if _has_active_srp_effector(args.dynamics_model.dynamic_effectors)
            throw(ArgumentError(
                "SimpleEphemeridesModel does not support solar-radiation-pressure ephemerides. " *
                "Use SpiceEphemeridesModel() for SRP-enabled runs."
            ))
        end
    end
    return nothing
end

function _initialize_heat_rate_buffers!(p)
    n_sats = length(p.args.dynamics_model.spacecraft)
    if length(p.shared_buffers.heat_rates) != n_sats
        resize!(p.shared_buffers.heat_rates, n_sats)
    end
    @inbounds for i in 1:n_sats
        n_links = length(p.args.dynamics_model.spacecraft[i].links)
        rates = p.shared_buffers.heat_rates[i]
        if !(rates isa Vector{Float64})
            rates = Float64[]
            p.shared_buffers.heat_rates[i] = rates
        end
        if length(rates) != n_links
            resize!(rates, n_links)
        end
        fill!(rates, 0.0)
    end
    return nothing
end

@inline function _gram_per_sat_instances_enabled()::Bool
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_GRAM_PER_SAT_INSTANCES", "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError(
        "Invalid SPACEAGORA_GRAM_PER_SAT_INSTANCES='$raw'. Use one of: 1/0, true/false, yes/no, on/off."
    ))
end

@inline function _parse_positive_float_env(name::String, default::Float64)::Float64
    raw = strip(_engine_env_get(name, string(default)))
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    parsed > 0.0 || throw(ArgumentError("$name must be > 0.0, got $parsed"))
    return parsed
end

@inline function _parse_unit_float_env(name::String, default::Float64)::Float64
    raw = strip(_engine_env_get(name, string(default)))
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    (0.0 < parsed <= 1.0) || throw(ArgumentError("$name must satisfy 0.0 < value <= 1.0, got $parsed"))
    return parsed
end

@inline function _parse_nonnegative_int_env(name::String, default::Int)::Int
    raw = strip(_engine_env_get(name, string(default)))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer value, got '$raw'"))
    end
    parsed >= 0 || throw(ArgumentError("$name must be >= 0, got $parsed"))
    return parsed
end

@inline function _srp_ephemeris_cache_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_SRP_EPHEMERIS_CACHE", true)
end

@inline function _srp_ephemeris_cache_dt_s()::Float64
    return _parse_positive_float_env("SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S", 30.0)
end

@inline function _srp_ephemeris_cache_max_samples()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES", 200_000)
end

@inline function _nbody_ephemeris_cache_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_NBODY_EPHEMERIS_CACHE", true)
end

@inline function _nbody_ephemeris_cache_dt_s()::Float64
    return _parse_positive_float_env("SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S", 30.0)
end

@inline function _nbody_ephemeris_cache_max_samples()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES", 200_000)
end

@inline function _planet_frame_cache_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_PLANET_FRAME_CACHE", true)
end

@inline function _planet_frame_cache_dt_s()::Float64
    return _parse_positive_float_env("SPACEAGORA_PLANET_FRAME_CACHE_DT_S", 30.0)
end

@inline function _planet_frame_cache_max_samples()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES", 400_000)
end

@inline function _spice_rhs_memo_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_SPICE_RHS_MEMO", true)
end

@inline function _ephemeris_reuse_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EPHEMERIS_CACHE_REUSE", true)
end

@inline function _ephemeris_reuse_max_entries()::Int
    return _parse_nonnegative_int_env("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES", 32)
end

@inline function _cache_time_key(x::Float64)::Int64
    return round(Int64, x * 1e6)
end

@inline function _planet_transform_key(planet)::NTuple{9, Int64}
    return ntuple(i -> round(Int64, planet.J2000_to_pci[i] * 1e12), 9)
end

@inline function _srp_ephemeris_reuse_key(primary_body_name::String, et_start::Float64, mission_end_s::Float64, dt_s::Float64)
    return (
        :srp,
        primary_body_name,
        _cache_time_key(et_start),
        _cache_time_key(mission_end_s),
        _cache_time_key(dt_s)
    )
end

@inline function _nbody_ephemeris_reuse_key(primary_body_name::String, body_query_names::Vector{String}, et_start::Float64, mission_end_s::Float64, dt_s::Float64)
    return (
        :nbody,
        primary_body_name,
        Tuple(body_query_names),
        _cache_time_key(et_start),
        _cache_time_key(mission_end_s),
        _cache_time_key(dt_s)
    )
end

@inline function _planet_frame_ephemeris_reuse_key(planet, ephemerides_model, et_start::Float64, mission_end_s::Float64, dt_s::Float64)
    return (
        :planet_frame,
        string(planet.name),
        SimulationModel.ephemerides_cache_key(ephemerides_model),
        _planet_transform_key(planet),
        _cache_time_key(et_start),
        _cache_time_key(mission_end_s),
        _cache_time_key(dt_s)
    )
end

@inline function _ephemeris_reuse_lookup(cache::Dict{Any, T}, key) where {T}
    return lock(_EPHEMERIS_REUSE_LOCK) do
        get(cache, key, nothing)
    end
end

@inline function _ephemeris_reuse_store!(cache::Dict{Any, T}, key, value::T, max_entries::Int)::T where {T}
    return lock(_EPHEMERIS_REUSE_LOCK) do
        existing = get(cache, key, nothing)
        if existing !== nothing
            return existing
        end
        if max_entries <= 0
            return value
        end
        if length(cache) >= max_entries
            oldest = first(keys(cache))
            delete!(cache, oldest)
        end
        cache[key] = value
        return value
    end
end

function _clear_ephemeris_reuse_cache!()
    lock(_EPHEMERIS_REUSE_LOCK) do
        empty!(_SRP_EPHEMERIS_REUSE_CACHE)
        empty!(_NBODY_EPHEMERIS_REUSE_CACHE)
        empty!(_PLANET_FRAME_EPHEMERIS_REUSE_CACHE)
    end
    return nothing
end

@inline function _has_active_srp_effector(dynamic_effectors::Tuple)::Bool
    @inbounds for effector in dynamic_effectors
        if effector isa SimulationModel.SolarRadiationPressureModel && effector.A > 0.0
            return true
        end
    end
    return false
end

@inline function _effector_parallel_mode()::Symbol
    return SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_EFFECTOR_PARALLEL"; default="auto")
end

@inline function _rhs_batch_parallel_mode()::Symbol
    return SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_RHS_BATCH_PARALLEL"; default="auto")
end

@inline function _profile_forces_serial_rhs()::Bool
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_PARALLEL_PROFILE", "")))
    return raw in ("r0", "serial", "r0_true_serial", "true_serial")
end

@inline function _rhs_batch_parallel_enabled(num_spacecraft::Int)::Bool
    if _profile_forces_serial_rhs()
        return false
    end
    mode = _rhs_batch_parallel_mode()
    if mode == :off
        return false
    elseif mode == :on
        return true
    end
    return num_spacecraft > 1 && Polyester.num_cores() > 1
end

@inline function _effector_thread_threshold()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_THREAD_THRESHOLD", 2)
end

@inline function _effector_max_threads()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_MAX_THREADS", 4)
end

@inline function _effector_outer_parallel_hint()::Bool
    return SimulationModel.ParallelPolicy.outer_parallel_active()
end

@inline function _effector_allow_with_outer()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EFFECTOR_PARALLEL_ALLOW_WITH_OUTER", false)
end

@inline function _effector_heavy_only()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY", true)
end

@inline function _effector_cost_ns_per_item_default()::Float64
    return _parse_positive_float_env("SPACEAGORA_EFFECTOR_COST_NS_PER_ITEM_DEFAULT", 2.5e4)
end

@inline function _effector_cost_min_samples()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_COST_MIN_SAMPLES", 4)
end

@inline function _effector_cost_ema_alpha()::Float64
    return _parse_unit_float_env("SPACEAGORA_EFFECTOR_COST_EMA_ALPHA", 0.2)
end

@inline function _effector_work_ns_per_worker_threshold()::Float64
    return _parse_positive_float_env("SPACEAGORA_EFFECTOR_WORK_NS_PER_WORKER_THRESHOLD", 4.0e4)
end

@inline function _effector_outer_work_scale()::Float64
    return _parse_positive_float_env("SPACEAGORA_EFFECTOR_OUTER_WORK_SCALE", 1.5)
end

@inline function _effector_long_mission_threshold_s()::Float64
    return _parse_positive_float_env("SPACEAGORA_EFFECTOR_LONG_MISSION_THRESHOLD_S", 5400.0)
end

@inline function _effector_long_orbit_threshold()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_LONG_ORBIT_THRESHOLD", 8)
end

@inline function _dynamic_effector_threadsafe(::Any)::Bool
    return false
end
@inline _dynamic_effector_threadsafe(::SimulationModel.InverseSquaredGravityModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.InverseSquaredJ2GravityModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.NBodyGravityModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.GravitationalHarmonicsModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.SolarRadiationPressureModel)::Bool = true
@inline _dynamic_effector_threadsafe(::SimulationModel.AerodynamicCoefficientfM)::Bool = true

@inline function _dynamic_effectors_parallel_supported(dynamic_effectors::Tuple)::Bool
    aero_fm_count = 0
    @inbounds for effector in dynamic_effectors
        if effector isa SimulationModel.AerodynamicCoefficientfM
            aero_fm_count += 1
        end
        _dynamic_effector_threadsafe(effector) || return false
    end
    return aero_fm_count <= 1
end

@inline function _mission_is_long_for_effector_threads(args)::Bool
    mission_cfg = args.mission_configuration
    if mission_cfg.mission_type == SimulationModel.MissionOrbits
        return mission_cfg.number_of_orbits >= _effector_long_orbit_threshold()
    end
    return mission_cfg.mission_time >= _effector_long_mission_threshold_s()
end

@inline function _effector_shared_buffers(p)
    if p === nothing || !hasproperty(p, :shared_buffers)
        return nothing
    end
    return getproperty(p, :shared_buffers)
end

@inline function _effector_observed_cost_ns_per_item(shared_buffers)::Float64
    default_cost = _effector_cost_ns_per_item_default()
    shared_buffers === nothing && return default_cost
    if !(hasproperty(shared_buffers, :effector_cost_ns_per_item) && hasproperty(shared_buffers, :effector_cost_samples))
        return default_cost
    end
    samples = Int(getproperty(shared_buffers, :effector_cost_samples)[])
    estimate = Float64(getproperty(shared_buffers, :effector_cost_ns_per_item)[])
    if samples >= _effector_cost_min_samples() && isfinite(estimate) && estimate > 0.0
        return estimate
    end
    return default_cost
end

@inline function _effector_satellite_share_budget(num_sats::Int, budget::Int)::Int
    sat_concurrency = max(1, min(max(1, num_sats), max(1, budget)))
    return max(1, fld(max(1, budget), sat_concurrency))
end

@inline function _update_effector_cost_model!(
    shared_buffers,
    n_effectors::Int,
    elapsed_ns::Int64,
    allotment::Int
)::Nothing
    shared_buffers === nothing && return nothing
    if !(hasproperty(shared_buffers, :effector_cost_ns_per_item) && hasproperty(shared_buffers, :effector_cost_samples))
        return nothing
    end
    n_effectors > 0 || return nothing
    wall_elapsed = max(1.0, Float64(elapsed_ns))
    # Scale by allotment to keep a rough "work per effector" estimate across serial/threaded samples.
    sample_ns_per_item = wall_elapsed * max(1, allotment) / max(1, n_effectors)
    α = _effector_cost_ema_alpha()
    estimate_ref = getproperty(shared_buffers, :effector_cost_ns_per_item)
    samples_ref = getproperty(shared_buffers, :effector_cost_samples)
    previous = Float64(estimate_ref[])
    estimate_ref[] = if isfinite(previous) && previous > 0.0
        (1.0 - α) * previous + α * sample_ns_per_item
    else
        sample_ns_per_item
    end
    samples_ref[] = min(typemax(Int64), samples_ref[] + Int64(1))
    return nothing
end

@inline function _dynamic_effector_thread_decision(
    args::SimulationConfiguration,
    p,
    dynamic_effectors::Tuple,
    num_sats::Int
)
    mode = _effector_parallel_mode()
    n_effectors = length(dynamic_effectors)
    if n_effectors <= 1
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    if !_dynamic_effectors_parallel_supported(dynamic_effectors)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end

    outer_active = _effector_outer_parallel_hint()
    allow_with_outer = _effector_allow_with_outer()
    budget = SimulationModel.ParallelPolicy.effective_inner_thread_budget()
    share_budget = _effector_satellite_share_budget(num_sats, budget)
    if outer_active && allow_with_outer
        share_budget = max(1, fld(share_budget, 2))
    end
    inner_floor = (!outer_active && num_sats > 1 && budget > 1) ? min(2, budget) : 1
    max_allotment = min(_effector_max_threads(), budget, max(share_budget, inner_floor))

    shared_buffers = _effector_shared_buffers(p)
    per_effector_cost_ns = _effector_observed_cost_ns_per_item(shared_buffers)
    estimated_work_ns = per_effector_cost_ns * n_effectors
    work_per_worker_ns = estimated_work_ns / max(1, max_allotment)
    target_ns = _effector_work_ns_per_worker_threshold() * (outer_active ? _effector_outer_work_scale() : 1.0)
    heavy_work = work_per_worker_ns >= target_ns

    policy = SimulationModel.ParallelPolicy.thread_policy_decision(
        n_effectors;
        mode=mode,
        threshold=_effector_thread_threshold(),
        heavy_work=heavy_work,
        heavy_only=_effector_heavy_only(),
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:dynamic_effectors
    )
    allotment = min(policy.allotment, max_allotment)
    use_threads = policy.use_threads && allotment > 1
    return (use_threads=use_threads, allotment=use_threads ? allotment : 1, mode=mode, policy_applied=true)
end

@inline function _dynamic_effector_thread_decision(
    args::SimulationConfiguration,
    dynamic_effectors::Tuple,
    num_sats::Int
)
    return _dynamic_effector_thread_decision(args, nothing, dynamic_effectors, num_sats)
end

function _initialize_density_model_instances!(p)
    instances = p.shared_buffers.density_models
    empty!(instances)

    density_model = p.args.environment_model.density_model
    if !_gram_per_sat_instances_enabled() || !(
        density_model isa SimulationModel.EnvironmentModels.GRAMAtmosphereModel ||
        density_model isa SimulationModel.EnvironmentModels.GRAMAtmosphereModelSurrogate
    )
        return nothing
    end

    n_sats = length(p.args.dynamics_model.spacecraft)
    sizehint!(instances, n_sats)
    @inbounds for _ in 1:n_sats
        # One GRAM handle per satellite avoids sharing mutable native model state.
        push!(instances, deepcopy(density_model))
    end
    return nothing
end

function _initialize_density_cache_buffers!(p)
    n_sats = length(p.args.dynamics_model.spacecraft)
    caches = p.shared_buffers.gram_density_cache
    if length(caches) != n_sats
        resize!(caches, n_sats)
    end
    fill!(caches, nothing)
    return nothing
end

function _initialize_gram_isolated_pool_buffers!(p)
    empty!(p.shared_buffers.gram_isolated_pool_models)
    empty!(p.shared_buffers.gram_isolated_pool_locks)
    return nothing
end

function _initialize_harmonics_workspace_buffers!(p)
    n_sats = length(p.args.dynamics_model.spacecraft)
    workspaces = p.shared_buffers.harmonics_workspaces
    if length(workspaces) != n_sats
        resize!(workspaces, n_sats)
    end
    fill!(workspaces, nothing)
    return nothing
end

function _initialize_nbody_workspace_buffers!(p)
    n_sats = length(p.args.dynamics_model.spacecraft)
    workspaces = p.shared_buffers.nbody_workspaces
    if length(workspaces) != n_sats
        resize!(workspaces, n_sats)
    end
    fill!(workspaces, nothing)
    return nothing
end

function _initialize_aero_workspace_buffers!(p)
    n_sats = length(p.args.dynamics_model.spacecraft)
    workspaces = p.shared_buffers.aero_workspaces
    if length(workspaces) != n_sats
        resize!(workspaces, n_sats)
    end
    fill!(workspaces, nothing)
    return nothing
end

function _initialize_srp_sun_cache_buffer!(p)
    p.shared_buffers.srp_sun_ephemeris_cache[] = nothing
    return nothing
end

function _initialize_nbody_ephemeris_cache_buffer!(p)
    p.shared_buffers.nbody_ephemeris_cache[] = nothing
    return nothing
end

function _initialize_planet_frame_cache_buffer!(p)
    p.shared_buffers.planet_frame_ephemeris_cache[] = nothing
    return nothing
end

function _initialize_spice_rhs_memo_mode!(p)
    p.shared_buffers.spice_rhs_memo_enabled[] = _spice_rhs_memo_enabled()
    return nothing
end

@inline function _reset_spice_runtime_counters!(p)
    counters = p.shared_buffers.spice_runtime_counters
    Base.Threads.atomic_xchg!(counters.nbody_spkpos_runtime_calls, 0)
    Base.Threads.atomic_xchg!(counters.nbody_spkpos_cache_build_calls, 0)
    Base.Threads.atomic_xchg!(counters.srp_spkpos_runtime_calls, 0)
    Base.Threads.atomic_xchg!(counters.srp_spkpos_cache_build_calls, 0)
    Base.Threads.atomic_xchg!(counters.planet_pxform_runtime_calls, 0)
    Base.Threads.atomic_xchg!(counters.planet_pxform_cache_build_calls, 0)
    return nothing
end

@inline function _reset_spice_rhs_memo!(p)
    memo = p.shared_buffers.spice_rhs_memo
    lock(memo.lock) do
        memo.et = NaN
        memo.primary_body_name = ""
        empty!(memo.body_positions_j2000)
    end
    return nothing
end

@inline function _spice_runtime_counters_snapshot(p)
    counters = p.shared_buffers.spice_runtime_counters
    nbody_runtime = counters.nbody_spkpos_runtime_calls[]
    nbody_build = counters.nbody_spkpos_cache_build_calls[]
    srp_runtime = counters.srp_spkpos_runtime_calls[]
    srp_build = counters.srp_spkpos_cache_build_calls[]
    planet_runtime = counters.planet_pxform_runtime_calls[]
    planet_build = counters.planet_pxform_cache_build_calls[]
    return (
        nbody_spkpos_runtime_calls=nbody_runtime,
        nbody_spkpos_cache_build_calls=nbody_build,
        nbody_spkpos_total_calls=(nbody_runtime + nbody_build),
        srp_spkpos_runtime_calls=srp_runtime,
        srp_spkpos_cache_build_calls=srp_build,
        srp_spkpos_total_calls=(srp_runtime + srp_build),
        planet_pxform_runtime_calls=planet_runtime,
        planet_pxform_cache_build_calls=planet_build,
        planet_pxform_total_calls=(planet_runtime + planet_build)
    )
end

@inline function _has_active_nbody_effector(dynamic_effectors::Tuple)::Bool
    @inbounds for effector in dynamic_effectors
        if effector isa SimulationModel.NBodyGravityModel && !isempty(effector.body_names)
            return true
        end
    end
    return false
end

function _collect_nbody_query_names(dynamic_effectors::Tuple)::Vector{String}
    names = String[]
    seen = Set{String}()
    @inbounds for effector in dynamic_effectors
        if !(effector isa SimulationModel.NBodyGravityModel)
            continue
        end
        for body_name in effector.body_names
            query_name = SimulationModel.DynamicEffectors._spice_query_name(body_name)
            if !(query_name in seen)
                push!(names, query_name)
                push!(seen, query_name)
            end
        end
    end
    return names
end

function _initialize_srp_sun_ephemeris_cache!(p, et_start::Float64, mission_end_s::Float64)
    _srp_ephemeris_cache_enabled() || return nothing
    _has_active_srp_effector(p.args.dynamics_model.dynamic_effectors) || return nothing
    if !isfinite(mission_end_s) || mission_end_s <= 0.0
        return nothing
    end

    dt_s = _srp_ephemeris_cache_dt_s()
    n_samples = max(2, Int(ceil(mission_end_s / dt_s)) + 1)
    max_samples = _srp_ephemeris_cache_max_samples()
    if n_samples > max_samples
        @warn "SRP ephemeris cache disabled: required samples=$n_samples exceeds SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES=$max_samples."
        return nothing
    end

    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(p.args.environment_model.planet.name)
    if _ephemeris_reuse_enabled()
        reuse_key = _srp_ephemeris_reuse_key(primary_body_name, et_start, mission_end_s, dt_s)
        reused = _ephemeris_reuse_lookup(_SRP_EPHEMERIS_REUSE_CACHE, reuse_key)
        if reused isa SimulationModel.SRPSunEphemerisCache
            p.shared_buffers.srp_sun_ephemeris_cache[] = reused
            return nothing
        end
    end
    ets = Vector{Float64}(undef, n_samples)
    positions = Vector{SVector{3, Float64}}(undef, n_samples)

    lock(SimulationModel.SPICE_LOCK) do
        @inbounds for sample_idx in 1:n_samples
            et = et_start + min((sample_idx - 1) * dt_s, mission_end_s)
            ets[sample_idx] = et
            positions[sample_idx] = SVector{3, Float64}(spkpos("sun", et, "J2000", "none", primary_body_name)[1])
            Base.Threads.atomic_add!(p.shared_buffers.spice_runtime_counters.srp_spkpos_cache_build_calls, 1)
        end
    end

    cache_value = SimulationModel.SRPSunEphemerisCache(ets, positions)
    if _ephemeris_reuse_enabled()
        reuse_key = _srp_ephemeris_reuse_key(primary_body_name, et_start, mission_end_s, dt_s)
        cache_value = _ephemeris_reuse_store!(_SRP_EPHEMERIS_REUSE_CACHE, reuse_key, cache_value, _ephemeris_reuse_max_entries())
    end
    p.shared_buffers.srp_sun_ephemeris_cache[] = cache_value
    if p.args.simulation_settings.verbose
        println(
            "Initialized SRP Sun ephemeris cache: samples=$(n_samples), dt=$(round(dt_s; digits=3)) s, span=$(round(mission_end_s; digits=3)) s."
        )
    end
    return nothing
end

function _initialize_nbody_ephemeris_cache!(p, et_start::Float64, mission_end_s::Float64)
    _nbody_ephemeris_cache_enabled() || return nothing
    _has_active_nbody_effector(p.args.dynamics_model.dynamic_effectors) || return nothing
    if !isfinite(mission_end_s) || mission_end_s <= 0.0
        return nothing
    end

    body_query_names = _collect_nbody_query_names(p.args.dynamics_model.dynamic_effectors)
    isempty(body_query_names) && return nothing

    dt_s = _nbody_ephemeris_cache_dt_s()
    n_samples = max(2, Int(ceil(mission_end_s / dt_s)) + 1)
    max_samples = _nbody_ephemeris_cache_max_samples()
    if n_samples > max_samples
        @warn "N-body ephemeris cache disabled: required samples=$n_samples exceeds SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES=$max_samples."
        return nothing
    end

    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(p.args.environment_model.planet.name)
    if _ephemeris_reuse_enabled()
        reuse_key = _nbody_ephemeris_reuse_key(primary_body_name, body_query_names, et_start, mission_end_s, dt_s)
        reused = _ephemeris_reuse_lookup(_NBODY_EPHEMERIS_REUSE_CACHE, reuse_key)
        if reused isa SimulationModel.NBodyEphemerisCache
            p.shared_buffers.nbody_ephemeris_cache[] = reused
            return nothing
        end
    end
    n_bodies = length(body_query_names)
    ets = Vector{Float64}(undef, n_samples)
    positions = Matrix{SVector{3, Float64}}(undef, n_samples, n_bodies)

    lock(SimulationModel.SPICE_LOCK) do
        @inbounds for sample_idx in 1:n_samples
            et = et_start + min((sample_idx - 1) * dt_s, mission_end_s)
            ets[sample_idx] = et
            for body_idx in 1:n_bodies
                body_query = body_query_names[body_idx]
                positions[sample_idx, body_idx] = SVector{3, Float64}(spkpos(body_query, et, "J2000", "none", primary_body_name)[1])
                Base.Threads.atomic_add!(p.shared_buffers.spice_runtime_counters.nbody_spkpos_cache_build_calls, 1)
            end
        end
    end

    body_index_by_name = Dict{String, Int}()
    @inbounds for (idx, body_name) in pairs(body_query_names)
        body_index_by_name[body_name] = idx
    end

    cache_value = SimulationModel.NBodyEphemerisCache(
        primary_body_name,
        body_query_names,
        body_index_by_name,
        ets,
        positions
    )
    if _ephemeris_reuse_enabled()
        reuse_key = _nbody_ephemeris_reuse_key(primary_body_name, body_query_names, et_start, mission_end_s, dt_s)
        cache_value = _ephemeris_reuse_store!(_NBODY_EPHEMERIS_REUSE_CACHE, reuse_key, cache_value, _ephemeris_reuse_max_entries())
    end
    p.shared_buffers.nbody_ephemeris_cache[] = cache_value
    if p.args.simulation_settings.verbose
        println(
            "Initialized N-body ephemeris cache: bodies=$(n_bodies), samples=$(n_samples), dt=$(round(dt_s; digits=3)) s, span=$(round(mission_end_s; digits=3)) s."
        )
    end
    return nothing
end

function _initialize_planet_frame_ephemeris_cache!(p, et_start::Float64, mission_end_s::Float64)
    _planet_frame_cache_enabled() || return nothing
    if !isfinite(mission_end_s) || mission_end_s <= 0.0
        return nothing
    end

    dt_s = _planet_frame_cache_dt_s()
    n_samples = max(2, Int(ceil(mission_end_s / dt_s)) + 1)
    max_samples = _planet_frame_cache_max_samples()
    if n_samples > max_samples
        @warn "Planet frame cache disabled: required samples=$n_samples exceeds SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES=$max_samples."
        return nothing
    end

    planet = p.args.environment_model.planet
    ephemerides_model = p.args.environment_model.ephemerides_model
    if _ephemeris_reuse_enabled()
        reuse_key = _planet_frame_ephemeris_reuse_key(planet, ephemerides_model, et_start, mission_end_s, dt_s)
        reused = _ephemeris_reuse_lookup(_PLANET_FRAME_EPHEMERIS_REUSE_CACHE, reuse_key)
        if reused isa SimulationModel.PlanetFrameEphemerisCache
            p.shared_buffers.planet_frame_ephemeris_cache[] = reused
            return nothing
        end
    end
    ets = Vector{Float64}(undef, n_samples)
    quaternions = Vector{SVector{4, Float64}}(undef, n_samples)

    @inbounds for sample_idx in 1:n_samples
        et = et_start + min((sample_idx - 1) * dt_s, mission_end_s)
        ets[sample_idx] = et
        l_pi = SimulationModel.planet_frame_lpi(planet, et, ephemerides_model)
        quaternions[sample_idx] = SimulationModel.dcm_to_quaternion(l_pi)
        if SimulationModel.ephemerides_requires_spice(ephemerides_model)
            Base.Threads.atomic_add!(p.shared_buffers.spice_runtime_counters.planet_pxform_cache_build_calls, 1)
        end
    end

    cache_value = SimulationModel.PlanetFrameEphemerisCache(ets, quaternions)
    if _ephemeris_reuse_enabled()
        reuse_key = _planet_frame_ephemeris_reuse_key(planet, ephemerides_model, et_start, mission_end_s, dt_s)
        cache_value = _ephemeris_reuse_store!(_PLANET_FRAME_EPHEMERIS_REUSE_CACHE, reuse_key, cache_value, _ephemeris_reuse_max_entries())
    end
    p.shared_buffers.planet_frame_ephemeris_cache[] = cache_value
    if p.args.simulation_settings.verbose
        println(
            "Initialized planet frame cache: samples=$(n_samples), dt=$(round(dt_s; digits=3)) s, span=$(round(mission_end_s; digits=3)) s."
        )
    end
    return nothing
end
