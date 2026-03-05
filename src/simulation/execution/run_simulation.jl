include(joinpath(@__DIR__, "..", "..", "utils", "Reference_system.jl"))
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

@inline _typed_normalize_warning_enabled() = get(ENV, "SPACEAGORA_WARN_NORMALIZE", "1") == "1"
@inline _typed_allow_legacy_normalize() = get(ENV, "SPACEAGORA_ALLOW_TYPED_NORMALIZE", "0") == "1"
@inline _typed_save_bundle_enabled() = get(ENV, "SPACEAGORA_SAVE_BUNDLE", "1") == "1"
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
    raw = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_PER_SAT_INSTANCES", "0")))
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
    raw = strip(get(ENV, name, string(default)))
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    parsed > 0.0 || throw(ArgumentError("$name must be > 0.0, got $parsed"))
    return parsed
end

@inline function _parse_unit_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    (0.0 < parsed <= 1.0) || throw(ArgumentError("$name must satisfy 0.0 < value <= 1.0, got $parsed"))
    return parsed
end

@inline function _parse_nonnegative_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
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

@inline function _planet_frame_ephemeris_reuse_key(planet, et_start::Float64, mission_end_s::Float64, dt_s::Float64)
    return (
        :planet_frame,
        string(planet.name),
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
    raw = lowercase(strip(get(ENV, "SPACEAGORA_PARALLEL_PROFILE", "")))
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
    if _ephemeris_reuse_enabled()
        reuse_key = _planet_frame_ephemeris_reuse_key(planet, et_start, mission_end_s, dt_s)
        reused = _ephemeris_reuse_lookup(_PLANET_FRAME_EPHEMERIS_REUSE_CACHE, reuse_key)
        if reused isa SimulationModel.PlanetFrameEphemerisCache
            p.shared_buffers.planet_frame_ephemeris_cache[] = reused
            return nothing
        end
    end
    frame_name = "IAU_$(planet.name)"
    ets = Vector{Float64}(undef, n_samples)
    quaternions = Vector{SVector{4, Float64}}(undef, n_samples)

    lock(SimulationModel.SPICE_LOCK) do
        @inbounds for sample_idx in 1:n_samples
            et = et_start + min((sample_idx - 1) * dt_s, mission_end_s)
            ets[sample_idx] = et
            l_pi = SMatrix{3, 3, Float64}(pxform("J2000", frame_name, et)) * planet.J2000_to_pci'
            quaternions[sample_idx] = SimulationModel.dcm_to_quaternion(l_pi)
            Base.Threads.atomic_add!(p.shared_buffers.spice_runtime_counters.planet_pxform_cache_build_calls, 1)
        end
    end

    cache_value = SimulationModel.PlanetFrameEphemerisCache(ets, quaternions)
    if _ephemeris_reuse_enabled()
        reuse_key = _planet_frame_ephemeris_reuse_key(planet, et_start, mission_end_s, dt_s)
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

@inline function _resolve_component_tolerance(component_tol::Float64, fallback_tol::Float64, name::String)::Float64
    if component_tol < 0.0
        throw(ArgumentError("$name must be >= 0.0, got $component_tol"))
    end
    return component_tol == 0.0 ? fallback_tol : component_tol
end

@inline function _requires_componentwise_tolerances(args)::Bool
    tol = args.integration_tolerances
    return args.mission_configuration.orientation_sim ||
           tol.reltol_mass != 0.0 || tol.abstol_mass != 0.0 ||
           tol.reltol_heat_load != 0.0 || tol.abstol_heat_load != 0.0
end

function _build_solver_tolerances(u_state::ComponentVector, args)
    tol = args.integration_tolerances
    if !_requires_componentwise_tolerances(args)
        return tol.reltol_orbit, tol.abstol_orbit
    end

    reltol_mass = _resolve_component_tolerance(tol.reltol_mass, tol.reltol_orbit, "reltol_mass")
    abstol_mass = _resolve_component_tolerance(tol.abstol_mass, tol.abstol_orbit, "abstol_mass")
    reltol_heat = _resolve_component_tolerance(tol.reltol_heat_load, tol.reltol_orbit, "reltol_heat_load")
    abstol_heat = _resolve_component_tolerance(tol.abstol_heat_load, tol.abstol_orbit, "abstol_heat_load")
    reltol_ω = _resolve_component_tolerance(tol.reltol_angular_rate, tol.reltol_orbit, "reltol_angular_rate")
    abstol_ω = _resolve_component_tolerance(tol.abstol_angular_rate, tol.abstol_orbit, "abstol_angular_rate")

    reltol_state = copy(u_state)
    abstol_state = copy(u_state)
    reltol_state .= tol.reltol_orbit
    abstol_state .= tol.abstol_orbit
    @inbounds for i in eachindex(reltol_state.sc)
        reltol_state.sc[i].mass = reltol_mass
        abstol_state.sc[i].mass = abstol_mass
        reltol_state.sc[i].heat_loads .= reltol_heat
        abstol_state.sc[i].heat_loads .= abstol_heat
    end

    if args.mission_configuration.orientation_sim
        @inbounds for i in eachindex(reltol_state.sc)
            reltol_state.sc[i].ω .= reltol_ω
            abstol_state.sc[i].ω .= abstol_ω
        end
    end

    if args.mission_configuration.orientation_sim
        @inbounds for i in eachindex(reltol_state.sc)
        reltol_state.sc[i].q .= tol.reltol_quaternion
        abstol_state.sc[i].q .= tol.abstol_quaternion
        end
    end
    return reltol_state, abstol_state
end

@inline function _solver_policy_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_SOLVER_MODE", "tsit5")))
    if mode in ("tsit5", "default")
        return :tsit5
    elseif mode in ("auto_stiff", "auto-stiff", "autostiff", "auto")
        return :auto_stiff
    elseif mode in ("rodas5p", "rodas", "stiff")
        return :rodas5p
    elseif mode in ("split_imex", "split-imex", "split", "imex")
        return :split_imex
    elseif mode in ("multirate", "multirate_split", "split_multirate", "mr")
        return :multirate
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_SOLVER_MODE='$mode'. Use one of: tsit5, auto_stiff, rodas5p, split_imex, multirate."
    ))
end

@inline function _retcode_is_stiff_symptom(retcode)::Bool
    rc = string(retcode)
    return rc in ("Unstable", "DtLessThanMin", "MaxIters", "InitialFailure")
end

@inline function _auto_stiff_switched(sol)::Bool
    hasproperty(sol, :alg_choice) || return false
    choices = getproperty(sol, :alg_choice)
    isempty(choices) && return false
    first_choice = first(choices)
    @inbounds for choice in choices
        if choice != first_choice
            return true
        end
    end
    return false
end

@inline function _solver_maxiters()::Union{Nothing, Int}
    raw = strip(get(ENV, "SPACEAGORA_SOLVER_MAXITERS", ""))
    isempty(raw) && return nothing
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be an integer, got '$raw'."))
    end
    parsed > 0 || throw(ArgumentError("SPACEAGORA_SOLVER_MAXITERS must be > 0, got $parsed."))
    return parsed
end

@inline function _split_imex_solver_spec()
    mode = lowercase(strip(get(ENV, "SPACEAGORA_SPLIT_IMEX_SOLVER", "kencarp4")))
    if mode in ("kencarp4", "ken4", "default")
        return (alg=KenCarp4(autodiff=AutoFiniteDiff()), label="KenCarp4")
    elseif mode in ("kencarp47", "ken47")
        return (alg=KenCarp47(autodiff=AutoFiniteDiff()), label="KenCarp47")
    elseif mode in ("kencarp58", "ken58")
        return (alg=KenCarp58(autodiff=AutoFiniteDiff()), label="KenCarp58")
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_SPLIT_IMEX_SOLVER='$mode'. Use one of: kencarp4, kencarp47, kencarp58."
    ))
end

@inline function _multirate_fast_substeps()::Int
    raw = strip(get(ENV, "SPACEAGORA_MULTIRATE_FAST_SUBSTEPS", "8"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be an integer, got '$raw'."))
    end
    parsed > 0 || throw(ArgumentError("SPACEAGORA_MULTIRATE_FAST_SUBSTEPS must be > 0, got $parsed."))
    return parsed
end

@inline function _multirate_slow_dt_s(args)::Float64
    default_dt = min(args.integration_tolerances.dt_max_orbit, 2.0)
    raw = strip(get(ENV, "SPACEAGORA_MULTIRATE_SLOW_DT_S", ""))
    dt = if isempty(raw)
        default_dt
    else
        parsed = try
            parse(Float64, raw)
        catch
            throw(ArgumentError("SPACEAGORA_MULTIRATE_SLOW_DT_S must be a number, got '$raw'."))
        end
        parsed
    end
    dt > 0.0 || throw(ArgumentError("SPACEAGORA_MULTIRATE_SLOW_DT_S must be > 0.0, got $dt."))
    return min(dt, args.integration_tolerances.dt_max_orbit)
end

@inline function _multirate_solver_spec(env_name::String, default_mode::String)
    mode = lowercase(strip(get(ENV, env_name, default_mode)))
    if mode in ("tsit5", "tsit", "default")
        return (alg=Tsit5(), label="Tsit5", auto_switch_capable=false)
    elseif mode in ("auto_stiff", "auto-stiff", "autostiff", "auto")
        return (
            alg=AutoTsit5(Rodas5P(autodiff=AutoFiniteDiff())),
            label="AutoTsit5(Rodas5P)",
            auto_switch_capable=true
        )
    elseif mode in ("rodas5p", "rodas", "stiff")
        return (alg=Rodas5P(autodiff=AutoFiniteDiff()), label="Rodas5P", auto_switch_capable=false)
    elseif mode in ("kencarp4", "ken4")
        return (alg=KenCarp4(autodiff=AutoFiniteDiff()), label="KenCarp4", auto_switch_capable=false)
    end
    throw(ArgumentError(
        "Unsupported $(env_name)='$mode'. Use one of: tsit5, auto_stiff, rodas5p, kencarp4."
    ))
end

@inline _multirate_slow_solver_spec() = _multirate_solver_spec("SPACEAGORA_MULTIRATE_SLOW_SOLVER", "tsit5")
@inline _multirate_fast_solver_spec() = _multirate_solver_spec("SPACEAGORA_MULTIRATE_FAST_SOLVER", "auto_stiff")

@inline function _solve_with_explicit_solver(prob, args, alg, reltol_tol, abstol_tol; dtmax_override::Union{Nothing, Float64}=nothing)
    maxiters = _solver_maxiters()
    dtmax_use = isnothing(dtmax_override) ? args.integration_tolerances.dt_max_orbit : dtmax_override
    dtmax_use > 0.0 || throw(ArgumentError("Solver dtmax must be > 0.0, got $dtmax_use."))
    if maxiters === nothing
        return solve(
            prob,
            alg;
            reltol=reltol_tol,
            abstol=abstol_tol,
            dtmax=dtmax_use
        )
    end
    return solve(
        prob,
        alg;
        reltol=reltol_tol,
        abstol=abstol_tol,
        dtmax=dtmax_use,
        maxiters=maxiters
    )
end

@inline function _split_subproblem(prob, f, u, tspan)
    return ODEProblem(f, u, tspan, prob.p; prob.kwargs...)
end

function _solve_with_multirate_solver(prob, args, reltol_tol, abstol_tol)
    if !(hasproperty(prob.f, :f1) && hasproperty(prob.f, :f2))
        throw(ArgumentError("SPACEAGORA_SOLVER_MODE=multirate requires a split problem with f1/f2 components."))
    end

    t_start = Float64(first(prob.tspan))
    t_end = Float64(last(prob.tspan))
    if t_end <= t_start
        sol = _solve_with_explicit_solver(prob, args, Tsit5(), reltol_tol, abstol_tol)
        return sol, (
            slow_solver="Tsit5",
            fast_solver="Tsit5",
            macro_steps=0,
            fast_substeps=0,
            slow_dt_s=0.0,
            fast_dt_s=0.0,
            auto_switch_events=0
        )
    end

    slow_spec = _multirate_slow_solver_spec()
    fast_spec = _multirate_fast_solver_spec()
    fast_substeps = _multirate_fast_substeps()
    slow_dt_s = _multirate_slow_dt_s(args)
    fast_dt_s = slow_dt_s / fast_substeps

    t_cursor = t_start
    u_cursor = deepcopy(prob.u0)
    final_sol = nothing
    macro_steps = 0
    auto_switch_events = 0

    while t_cursor < t_end
        t_next = min(t_cursor + slow_dt_s, t_end)
        macro_steps += 1

        # Strang splitting: fast half-step -> slow full-step -> fast half-step.
        segment_dt = t_next - t_cursor
        half_dt = 0.5 * segment_dt
        t_half = t_cursor + half_dt

        if half_dt > 0.0
            fast_prob_pre = _split_subproblem(prob, prob.f.f2, u_cursor, (t_cursor, t_half))
            sol_fast_pre = _solve_with_explicit_solver(
                fast_prob_pre,
                args,
                fast_spec.alg,
                reltol_tol,
                abstol_tol;
                dtmax_override=min(fast_dt_s, half_dt)
            )
            if fast_spec.auto_switch_capable && _auto_stiff_switched(sol_fast_pre)
                auto_switch_events += 1
            end
            if !SciMLBase.successful_retcode(sol_fast_pre.retcode)
                return sol_fast_pre, (
                    slow_solver=slow_spec.label,
                    fast_solver=fast_spec.label,
                    macro_steps=macro_steps,
                    fast_substeps=fast_substeps,
                    slow_dt_s=slow_dt_s,
                    fast_dt_s=fast_dt_s,
                    auto_switch_events=auto_switch_events
                )
            end
            u_cursor = deepcopy(sol_fast_pre.u[end])
            final_sol = sol_fast_pre
        end

        slow_prob = _split_subproblem(prob, prob.f.f1, u_cursor, (t_cursor, t_next))
        sol_slow = _solve_with_explicit_solver(
            slow_prob,
            args,
            slow_spec.alg,
            reltol_tol,
            abstol_tol;
            dtmax_override=segment_dt
        )
        if slow_spec.auto_switch_capable && _auto_stiff_switched(sol_slow)
            auto_switch_events += 1
        end
        if !SciMLBase.successful_retcode(sol_slow.retcode)
            return sol_slow, (
                slow_solver=slow_spec.label,
                fast_solver=fast_spec.label,
                macro_steps=macro_steps,
                fast_substeps=fast_substeps,
                slow_dt_s=slow_dt_s,
                fast_dt_s=fast_dt_s,
                auto_switch_events=auto_switch_events
            )
        end
        u_cursor = deepcopy(sol_slow.u[end])
        final_sol = sol_slow

        if half_dt > 0.0
            fast_prob_post = _split_subproblem(prob, prob.f.f2, u_cursor, (t_half, t_next))
            sol_fast_post = _solve_with_explicit_solver(
                fast_prob_post,
                args,
                fast_spec.alg,
                reltol_tol,
                abstol_tol;
                dtmax_override=min(fast_dt_s, half_dt)
            )
            if fast_spec.auto_switch_capable && _auto_stiff_switched(sol_fast_post)
                auto_switch_events += 1
            end
            if !SciMLBase.successful_retcode(sol_fast_post.retcode)
                return sol_fast_post, (
                    slow_solver=slow_spec.label,
                    fast_solver=fast_spec.label,
                    macro_steps=macro_steps,
                    fast_substeps=fast_substeps,
                    slow_dt_s=slow_dt_s,
                    fast_dt_s=fast_dt_s,
                    auto_switch_events=auto_switch_events
                )
            end
            u_cursor = deepcopy(sol_fast_post.u[end])
            final_sol = sol_fast_post
        end

        t_cursor = t_next
    end

    return final_sol, (
        slow_solver=slow_spec.label,
        fast_solver=fast_spec.label,
        macro_steps=macro_steps,
        fast_substeps=fast_substeps,
        slow_dt_s=slow_dt_s,
        fast_dt_s=fast_dt_s,
        auto_switch_events=auto_switch_events
    )
end

function _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
    mode = _solver_policy_mode()
    if mode == :rodas5p
        sol = _solve_with_explicit_solver(prob, args, Rodas5P(autodiff=AutoFiniteDiff()), reltol_tol, abstol_tol)
        return sol, (
            solver="Rodas5P",
            initial_solver="Rodas5P",
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :auto_stiff
        # True stiffness-aware autoswitching handled internally by OrdinaryDiffEq.
        # This replaces the manual "retry with Rodas5P on Tsit5 failure" policy.
        autoswitch_alg = AutoTsit5(Rodas5P(autodiff=AutoFiniteDiff()))
        sol = _solve_with_explicit_solver(prob, args, autoswitch_alg, reltol_tol, abstol_tol)
        switched = _auto_stiff_switched(sol)
        return sol, (
            solver="AutoTsit5(Rodas5P)",
            initial_solver="AutoTsit5",
            fallback_used=switched,
            trigger_retcode=switched ? "internal_autoswitch" : missing
        )
    end

    if mode == :split_imex
        split_solver = _split_imex_solver_spec()
        sol = _solve_with_explicit_solver(prob, args, split_solver.alg, reltol_tol, abstol_tol)
        return sol, (
            solver="$(split_solver.label)(IMEX)",
            initial_solver=split_solver.label,
            fallback_used=false,
            trigger_retcode=missing
        )
    end

    if mode == :multirate
        sol, multirate_meta = _solve_with_multirate_solver(prob, args, reltol_tol, abstol_tol)
        switched = multirate_meta.auto_switch_events > 0
        return sol, (
            solver="Multirate(Strang; slow=$(multirate_meta.slow_solver), fast=$(multirate_meta.fast_solver))",
            initial_solver=multirate_meta.slow_solver,
            fallback_used=switched,
            trigger_retcode=switched ? "internal_autoswitch" : missing
        )
    end

    tsit_sol = _solve_with_explicit_solver(prob, args, Tsit5(), reltol_tol, abstol_tol)
    return tsit_sol, (
        solver="Tsit5",
        initial_solver="Tsit5",
        fallback_used=false,
        trigger_retcode=missing
    )
end

function _warn_legacy_normalize_flag!(args)
    if !args.simulation_settings.normalize || !_typed_normalize_warning_enabled()
        return nothing
    end
    if _normalize_warning_emitted[]
        return nothing
    end
    _normalize_warning_emitted[] = true
    @warn "SimulationSettings.normalize=true is legacy-only in typed run_simulation; propagation is always SI-native (m, s, kg). Set normalize=false to silence this warning."
    return nothing
end

function _enforce_typed_normalize_policy!(args)
    if !args.simulation_settings.normalize
        return nothing
    end
    if _typed_allow_legacy_normalize()
        _warn_legacy_normalize_flag!(args)
        return nothing
    end
    throw(ArgumentError(
        "SimulationSettings.normalize=true is unsupported in typed run_simulation. " *
        "Typed propagation is SI-native (m, s, kg). Set normalize=false, or set " *
        "SPACEAGORA_ALLOW_TYPED_NORMALIZE=1 only for legacy transition checks."
    ))
end

function _debug_print_nan_parameter_paths!(x, path::AbstractString="p")
    if x isa Number
        if isnan(x)
            println("NaN found in parameter: $path")
        end
        return nothing
    end

    if x isa Base.RefValue{<:Number}
        xv = x[]
        if isnan(xv)
            println("NaN found in parameter: $path[]")
        end
        return nothing
    end

    if x isa AbstractArray{<:Number}
        for (idx, xv) in pairs(x)
            if isnan(xv)
                println("NaN found in parameter: $path[$idx]")
            end
        end
        return nothing
    end

    # Skip generic arrays of non-numeric types to keep debug scans bounded.
    if x isa AbstractArray
        return nothing
    end

    T = typeof(x)
    if isstructtype(T)
        for field in fieldnames(T)
            val = getfield(x, field)
            _debug_print_nan_parameter_paths!(val, string(path, ".", field))
        end
    end
    return nothing
end

@inline function _results_bundle_prefix(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results")
end

@inline function _legacy_results_csv_path(args)::String
    return joinpath(args.simulation_settings.results_directory, "simulation_results.csv")
end

@inline function _collision_results_csv_path(args)::String
    stamp = Dates.format(now(UTC), dateformat"yyyymmddTHHMMSSsss")
    token = string(stamp, ".", getpid(), ".", rand(UInt))
    return joinpath(args.simulation_settings.results_directory, "simulation_results.$token.csv")
end

function _atomic_write_file(path::String, writer::Function; force::Bool=true)
    dir = dirname(path)
    mkpath(dir)
    tmp = joinpath(
        dir,
        string(".", basename(path), ".tmp.", getpid(), ".", Threads.threadid(), ".", rand(UInt))
    )
    try
        writer(tmp)
        mv(tmp, path; force=force)
    finally
        if isfile(tmp)
            rm(tmp; force=true)
        end
    end
    return path
end

function _write_legacy_results_csv!(results_df::DataFrame, args)::String
    primary_path = _legacy_results_csv_path(args)
    started_s = time()
    existed_before = isfile(primary_path)
    try
        return _atomic_write_file(primary_path, tmp -> CSV.write(tmp, results_df); force=false)
    catch err
        if err isa ArgumentError && isfile(primary_path)
            mtime_s = try
                stat(primary_path).mtime
            catch
                0.0
            end
            # Keep per-run collision artifacts only when a concurrent writer is likely.
            concurrent_writer = (!existed_before) || (mtime_s >= started_s)
            if concurrent_writer
                collision_path = _collision_results_csv_path(args)
                _atomic_write_file(collision_path, tmp -> CSV.write(tmp, results_df); force=false)
            end
            return _atomic_write_file(primary_path, tmp -> CSV.write(tmp, results_df); force=true)
        end
        rethrow(err)
    end
end

function _sha256_hex(path::String)::String
    open(path, "r") do io
        return bytes2hex(SHA.sha256(read(io)))
    end
end

@inline function _checkpoint_directory(args)::String
    if isempty(args.simulation_settings.checkpoint_directory)
        return joinpath(args.simulation_settings.results_directory, "checkpoints")
    end
    return args.simulation_settings.checkpoint_directory
end

@inline function _checkpoint_paths(args)
    ckpt_dir = _checkpoint_directory(args)
    return (
        data=joinpath(ckpt_dir, "simulation_checkpoint.bin"),
        manifest=joinpath(ckpt_dir, "simulation_checkpoint.manifest.toml")
    )
end

function _write_checkpoint!(args, t::Float64, u_state)
    paths = _checkpoint_paths(args)
    payload = (
        schema_version=CHECKPOINT_SCHEMA_VERSION,
        created_utc=string(now(UTC)),
        t=t,
        u=deepcopy(u_state)
    )
    _atomic_write_file(paths.data, tmp -> open(tmp, "w") do io
        serialize(io, payload)
    end)

    manifest = Dict{String, Any}(
        "schema_version" => CHECKPOINT_SCHEMA_VERSION,
        "created_utc" => string(now(UTC)),
        "time_s" => t,
        "data_path" => paths.data,
        "data_size_bytes" => filesize(paths.data),
        "data_sha256" => _sha256_hex(paths.data)
    )
    _atomic_write_file(paths.manifest, tmp -> open(tmp, "w") do io
        TOML.print(io, manifest)
    end)
    return nothing
end

function _load_checkpoint(args)
    paths = _checkpoint_paths(args)
    if !isfile(paths.data)
        return nothing
    end
    payload = open(paths.data, "r") do io
        deserialize(io)
    end
    if !haskey(payload, :t) || !haskey(payload, :u)
        throw(ArgumentError("Checkpoint payload missing required keys (:t, :u)."))
    end
    return (t=Float64(payload[:t]), u=payload[:u], data_path=paths.data, manifest_path=paths.manifest)
end

function _clear_checkpoint!(args)
    paths = _checkpoint_paths(args)
    isfile(paths.data) && rm(paths.data; force=true)
    isfile(paths.manifest) && rm(paths.manifest; force=true)
    return nothing
end

function _append_saved_segment!(times_acc::Vector{Float64}, data_acc::Vector, saved_values)
    seg_len = length(saved_values.t)
    seg_len == 0 && return nothing
    start_idx = 1
    if !isempty(times_acc) && isapprox(times_acc[end], saved_values.t[1]; atol=0.0, rtol=0.0)
        start_idx = 2
    end
    if start_idx <= seg_len
        append!(times_acc, @view saved_values.t[start_idx:seg_len])
        append!(data_acc, @view saved_values.saveval[start_idx:seg_len])
    end
    return nothing
end

@inline function _is_flat_scalar(value)::Bool
    return value === missing || value === nothing || value isa Number || value isa AbstractString || value isa Symbol || value isa Bool
end

function _find_sample_value(series)
    for value in series
        if value !== nothing
            return value
        end
    end
    return nothing
end

function _append_series_columns!(results_df::DataFrame, prefix::String, series)
    sample = _find_sample_value(series)
    if sample === nothing || _is_flat_scalar(sample)
        results_df[!, prefix] = collect(series)
        return nothing
    end

    if sample isa NamedTuple
        for key in keys(sample)
            child_series = [value === nothing ? nothing : getproperty(value, key) for value in series]
            _append_series_columns!(results_df, string(prefix, "_", key), child_series)
        end
        return nothing
    end

    if sample isa AbstractDict
        for key in sort!(collect(keys(sample)); by=string)
            child_series = [value === nothing ? nothing : value[key] for value in series]
            _append_series_columns!(results_df, string(prefix, "_", key), child_series)
        end
        return nothing
    end

    if sample isa Tuple || sample isa AbstractArray
        for idx in eachindex(sample)
            child_series = [value === nothing ? nothing : value[idx] for value in series]
            _append_series_columns!(results_df, string(prefix, "_", idx), child_series)
        end
        return nothing
    end

    results_df[!, prefix] = collect(series)
    return nothing
end

function _append_save_field_columns!(results_df::DataFrame, field, saved_data::Vector, num_sats::Int)
    field_series = [snapshot[field.name] for snapshot in saved_data]
    if field.per_satellite
        for sat_idx in 1:num_sats
            sat_series = [value[sat_idx] for value in field_series]
            _append_series_columns!(results_df, "sc$(sat_idx)_$(field.column_prefix)", sat_series)
        end
        return nothing
    end
    _append_series_columns!(results_df, field.column_prefix, field_series)
    return nothing
end

function _build_results_dataframe(times::Vector{Float64}, saved_data::Vector, save_fields, args)::DataFrame
    results_df = DataFrame(time=times)
    num_sats = length(args.dynamics_model.spacecraft)
    for field in save_fields
        _append_save_field_columns!(results_df, field, saved_data, num_sats)
    end
    return results_df
end

function _write_results_bundle!(
    results_df::DataFrame,
    times::Vector{Float64},
    args;
    csv_path::Union{Nothing, String}=nothing
)
    prefix = _results_bundle_prefix(args)
    feather_path = prefix * ".feather"
    manifest_path = prefix * ".manifest.toml"

    _atomic_write_file(feather_path, tmp -> Arrow.write(tmp, results_df))

    files = Dict{String, Any}()
    files["feather"] = Dict(
        "path" => feather_path,
        "size_bytes" => filesize(feather_path),
        "sha256" => _sha256_hex(feather_path)
    )

    if args.simulation_settings.save_csv
        csv_file_path = csv_path === nothing ? (prefix * ".csv") : csv_path
        if csv_path === nothing
            _atomic_write_file(csv_file_path, tmp -> CSV.write(tmp, results_df))
        end
        files["csv"] = Dict(
            "path" => csv_file_path,
            "size_bytes" => filesize(csv_file_path),
            "sha256" => _sha256_hex(csv_file_path)
        )
    end

    manifest = Dict{String, Any}(
        "schema_version" => RESULTS_BUNDLE_SCHEMA_VERSION,
        "created_utc" => string(now(UTC)),
        "mission_time_s" => args.mission_configuration.mission_time,
        "steps" => length(times),
        "spacecraft_count" => length(args.dynamics_model.spacecraft),
        "orientation_sim" => args.mission_configuration.orientation_sim,
        "files" => files
    )

    _atomic_write_file(manifest_path, tmp -> begin
        open(tmp, "w") do io
            TOML.print(io, manifest)
        end
    end)

    return nothing
end

@inline function _build_typed_solver_problem(u0, tspan, p, callbacks)
    mode = _solver_policy_mode()
    if mode == :split_imex || mode == :multirate
        return SplitODEProblem(
            spacecraft_dynamics_slow!,
            spacecraft_dynamics_fast_control!,
            u0,
            tspan,
            p,
            callback=callbacks
        )
    end
    return ODEProblem(spacecraft_dynamics!, u0, tspan, p, callback=callbacks)
end

function run_simulation(
    args;
    isolate_state::Bool=true,
    return_solution::Bool=false,
    return_solver_metadata::Bool=false,
    save_fields=nothing
)
    return SimulationModel.ParallelPolicy.with_policy_context() do
    # Isolate mutable campaign/model state by default so repeated/concurrent runs
    # do not alias shared in-memory objects.
    args = isolate_state ? deepcopy(args) : args

    # Typed pipeline is SI-native (meters, seconds, kilograms). The
    # `simulation_settings.normalize` field is legacy-only and rejected by default.
    _enforce_typed_normalize_policy!(args)
    _validate_orientation_inertia!(args)
    _validate_thermal_model_support!(args)
    try
        SimulationModel.ParallelPolicy.reset_policy_telemetry!()
        if SimulationModel.ParallelPolicy.persistent_hints_state_reset_requested()
            SimulationModel.ParallelPolicy.reset_persistent_hint_state!()
        end
    catch
    end

    # Set up the model and initial conditions
    initial_conditions = build_initial_conditions(args)
    if args.simulation_settings.verbose
        println("Initial conditions:")
        println(initial_conditions)
    end

    # Define the ODE parameters and callbacks
    p = SimulationModel.ODEParams{length(args.dynamics_model.spacecraft)}(args=args) # Define the parameters for the ODE problem, including the shared buffers for the callbacks
    _initialize_heat_rate_buffers!(p)
    _initialize_density_model_instances!(p)
    _initialize_density_cache_buffers!(p)
    _initialize_gram_isolated_pool_buffers!(p)
    _initialize_harmonics_workspace_buffers!(p)
    _initialize_nbody_workspace_buffers!(p)
    _initialize_aero_workspace_buffers!(p)
    _initialize_nbody_ephemeris_cache_buffer!(p)
    _initialize_srp_sun_cache_buffer!(p)
    _initialize_planet_frame_cache_buffer!(p)
    _initialize_spice_rhs_memo_mode!(p)
    _reset_spice_runtime_counters!(p)
    _reset_spice_rhs_memo!(p)
    p.shared_buffers.debug_control[] = get(ENV, "SPACEAGORA_DEBUG_CONTROL", "0") == "1"
    p.shared_buffers.debug_initial_derivative[] = get(ENV, "SPACEAGORA_DEBUG_INITIAL_DERIVATIVE", "0") == "1"
    save_fields_resolved = isnothing(save_fields) ? SimulationModel.default_save_fields(args) : collect(save_fields)
    save_field_names = Symbol[field.name for field in save_fields_resolved]
    length(unique(save_field_names)) == length(save_field_names) || throw(ArgumentError("save_fields names must be unique. Got $(save_field_names)."))
    saved_values = SavedValues(Float64, SimulationModel.SaveData)
    callbacks = SimulationModel.get_callbacks(
        length(args.dynamics_model.spacecraft),
        args.dynamics_model.dynamic_effectors,
        args;
        saved_values=saved_values,
        save_fields=save_fields_resolved
    ) # Get the callbacks based on the number of satellites and the dynamic effectors being used in the simulation
    initial_time = args.initial_time
    start_epoch = from_utc(DateTime(
            initial_time.year,
            initial_time.month,
            initial_time.day,
            initial_time.hour,
            initial_time.minute,
            initial_time.second
        ))
    et_start = lock(SimulationModel.SPICE_LOCK) do
        utc2et(to_utc(start_epoch))
    end
    p.shared_buffers.et_start[] = et_start
    lock(SimulationModel.SPICE_LOCK) do
        args.environment_model.planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(args.environment_model.planet.name)", et_start)) * args.environment_model.planet.J2000_to_pci' # Initialize the planet frame at the start of the simulation (will be updated in the callback)
    end
    Base.Threads.atomic_add!(p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls, 1)
    mission_end = args.mission_configuration.mission_time
    _initialize_nbody_ephemeris_cache!(p, et_start, mission_end)
    _initialize_srp_sun_ephemeris_cache!(p, et_start, mission_end)
    _initialize_planet_frame_ephemeris_cache!(p, et_start, mission_end)
    checkpoint_active = _typed_checkpoint_enabled(args)
    if checkpoint_active && args.simulation_settings.checkpoint_interval_s <= 0.0
        throw(ArgumentError("SimulationSettings.checkpoint_interval_s must be > 0 when checkpointing is enabled."))
    end

    t_start = 0.0
    u_start = initial_conditions
    if args.simulation_settings.resume_from_checkpoint
        ckpt = _load_checkpoint(args)
        if ckpt === nothing
            @warn "resume_from_checkpoint=true but no checkpoint file was found; starting from initial conditions."
        else
            t_start = ckpt.t
            u_start = ckpt.u
            if args.simulation_settings.verbose
                println("Resuming simulation from checkpoint at t=$(round(t_start, digits=6)) s")
            end
        end
    end

    # println("Initial conditions:")
    # println(initial_conditions)
    # println("ODE parameters:")
    # println(p)
    # println("args.mission_configuration.mission_time: $(args.mission_configuration.mission_time)")
    prob_debug = ODEProblem(spacecraft_dynamics!, u_start, (t_start, mission_end), p, callback=callbacks)
    if p.shared_buffers.debug_initial_derivative[]
        # 1. Manually evaluate the derivative at the start
        du_test = copy(prob_debug.u0)
        try
            prob_debug.f(du_test, prob_debug.u0, prob_debug.p, prob_debug.tspan[1])
        catch e
            @error "The derivative function itself crashed!" exception=e
            throw(ErrorException("Initial derivative evaluation failed; aborting solve in debug mode."))
        end

        # 2. Check for NaNs and print exactly where they are
        if any(isnan, du_test)
            println("--- INITIAL NaN DETECTED ---")

            # Check global parameters in p
            _debug_print_nan_parameter_paths!(prob_debug.p)

            # Check the state vector (u)
            # Assuming your u has a .sc field for satellites
            for (i, sat) in enumerate(du_test.sc)
                if any(isnan, sat.pos) || any(isnan, sat.vel)
                    println("NaN found in Satellite $i derivative!")
                    println("  Pos: $(sat.pos)")
                    println("  Vel: $(sat.vel)")
                end
            end

            throw(ErrorException("Initial derivative contains NaN; aborting solve in debug mode."))
        end
    end

    reltol_tol, abstol_tol = _build_solver_tolerances(u_start, args)
    last_sol = nothing
    solver_trace = NamedTuple[]
    checkpoint_saved_times = Float64[]
    checkpoint_saved_data = SimulationModel.SaveData[]

    if t_start < mission_end && checkpoint_active
        interval = args.simulation_settings.checkpoint_interval_s
        t_cursor = t_start
        u_cursor = deepcopy(u_start)

        while t_cursor < mission_end
            t_next = min(t_cursor + interval, mission_end)
            empty!(saved_values.t)
            empty!(saved_values.saveval)
            prob = _build_typed_solver_problem(u_cursor, (t_cursor, t_next), p, callbacks)
            seg_sol, solve_meta = _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
            push!(solver_trace, solve_meta)
            if !SciMLBase.successful_retcode(seg_sol.retcode)
                throw(ErrorException("Checkpointed solve failed with retcode=$(seg_sol.retcode)."))
            end
            _append_saved_segment!(checkpoint_saved_times, checkpoint_saved_data, saved_values)
            last_sol = seg_sol
            t_cursor = Float64(seg_sol.t[end])
            u_cursor = deepcopy(seg_sol.u[end])
            _write_checkpoint!(args, t_cursor, u_cursor)
        end

    elseif t_start < mission_end
        prob = _build_typed_solver_problem(u_start, (t_start, mission_end), p, callbacks)
        sol, solve_meta = _solve_with_solver_policy(prob, args, reltol_tol, abstol_tol)
        push!(solver_trace, solve_meta)
        if !SciMLBase.successful_retcode(sol.retcode)
            throw(ErrorException("Solve failed with retcode=$(sol.retcode)."))
        end
        last_sol = sol
    end

    # Process and save results
    if args.simulation_settings.results
        results_times = checkpoint_active ? checkpoint_saved_times : saved_values.t
        results_data = checkpoint_active ? checkpoint_saved_data : saved_values.saveval
        results_df = _build_results_dataframe(results_times, results_data, save_fields_resolved, args)
        # Keep backwards-compatible CSV contract used by existing scripts/tests.
        csv_path = _write_legacy_results_csv!(results_df, args)
        if _typed_save_bundle_enabled()
            _write_results_bundle!(results_df, results_times, args; csv_path=csv_path)
        end
    end

    if return_solution
        if checkpoint_active && args.simulation_settings.checkpoint_interval_s < mission_end
            @warn "return_solution=true with checkpointed integration returns the final segment ODESolution, not a stitched full-history ODESolution."
        end
        if return_solver_metadata
            parallel_policy = try
                SimulationModel.ParallelPolicy.policy_telemetry_snapshot()
            catch
                nothing
            end
            return (
                solution=last_sol,
                solver_mode=string(_solver_policy_mode()),
                solver_trace=solver_trace,
                parallel_policy=parallel_policy,
                spice_counters=_spice_runtime_counters_snapshot(p)
            )
        end
        return last_sol
    end
    return nothing
    end
end

@inline function _accumulate_dynamic_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    dynamic_effectors::Tuple,
    effector_decision
)
    effector_started_ns = time_ns()
    n_effectors = length(dynamic_effectors)
    if effector_decision.use_threads
        reduced = SimulationModel.ParallelPolicy.threaded_reduce(
            n_effectors,
            effector_decision.allotment,
            () -> MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (local_sum, eff_idx) -> begin
                effector = dynamic_effectors[eff_idx]
                force, torque = SimulationModel.calcForceTorque(effector, sc_view, p, sat_idx)
                local_sum[1] += force[1]
                local_sum[2] += force[2]
                local_sum[3] += force[3]
                local_sum[4] += torque[1]
                local_sum[5] += torque[2]
                local_sum[6] += torque[3]
                return nothing
            end,
            (dest, src) -> begin
                @inbounds for i in 1:6
                    dest[i] += src[i]
                end
                return nothing
            end
        )
        forces .= SVector{3, Float64}(reduced[1], reduced[2], reduced[3])
        torques .= SVector{3, Float64}(reduced[4], reduced[5], reduced[6])
    else
        @inbounds for effector in dynamic_effectors
            force, torque = SimulationModel.calcForceTorque(effector, sc_view, p, sat_idx)
            forces .+= force
            torques .+= torque
        end
    end
    elapsed_ns = Int64(time_ns() - effector_started_ns)
    if effector_decision.policy_applied && sat_idx == 1
        _update_effector_cost_model!(p.shared_buffers, n_effectors, elapsed_ns, effector_decision.allotment)
    end
    if effector_decision.policy_applied
        SimulationModel.ParallelPolicy.record_policy_observation!(
            :dynamic_effectors;
            mode=effector_decision.mode,
            num_items=n_effectors,
            use_threads=effector_decision.use_threads,
            elapsed_ns=elapsed_ns
        )
    end
    return nothing
end

@inline function _accumulate_control_effectors!(
    forces::MVector{3, Float64},
    torques::MVector{3, Float64},
    sc_view,
    p,
    sat_idx::Int,
    t::Float64,
    debug_control::Bool
)::Float64
    mass_rate = 0.0
    @inbounds for control_effector in p.args.control_model.control_effectors
        control_force, control_torque = SimulationModel.calcControlForceTorque(control_effector, sc_view, p, sat_idx, t)
        control_mass_rate = SimulationModel.calcControlMassFlowRate(control_effector, sc_view, p, sat_idx, t)
        if debug_control && (norm(control_force) > 0.0 || norm(control_torque) > 0.0)
            println("Applying control effect for spacecraft $sat_idx at time $t seconds:")
            println("  Control force: $control_force")
        end
        forces .+= control_force
        torques .+= control_torque
        mass_rate += isfinite(control_mass_rate) ? control_mass_rate : 0.0
    end
    return mass_rate
end

@inline function _assign_heat_rate_derivative!(du_heat::AbstractVector, heat_rates::AbstractVector)
    if length(heat_rates) == length(du_heat)
        du_heat .= heat_rates
        return nothing
    end
    du_heat .= 0.0
    n_copy = min(length(heat_rates), length(du_heat))
    @inbounds for j in 1:n_copy
        du_heat[j] = heat_rates[j]
    end
    return nothing
end

function spacecraft_dynamics!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    effector_decision = _dynamic_effector_thread_decision(p.args, p, dynamic_effectors, length(spacecraft))
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
            end
        end
    end
end # function spacecraft_dynamics!

function spacecraft_dynamics_slow!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    dynamics_model = p.args.dynamics_model
    dynamic_effectors = dynamics_model.dynamic_effectors
    spacecraft = dynamics_model.spacecraft
    p.shared_buffers.current_time[] = t
    effector_decision = _dynamic_effector_thread_decision(p.args, p, dynamic_effectors, length(spacecraft))
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = 0.0

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                _accumulate_dynamic_effectors!(forces, torques, sc_view, p, i, dynamic_effectors, effector_decision)

                du_view.pos .= sc_view.vel
                du_view.vel .= forces / sc_view.mass
                du_view.mass = 0.0

                if p.args.mission_configuration.orientation_sim
                    ω_body = SVector{3, Float64}(sc_view.ω)
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.5 * quat_mult(SVector{4, Float64}(ω_body..., 0.0), sc_view.q)
                    du_view.ω .= inertia_tensor \ (τ_body - cross(ω_body, inertia_tensor * ω_body))
                end

                _assign_heat_rate_derivative!(du_view.heat_loads, p.shared_buffers.heat_rates[i])
            end
        end
    end
end # function spacecraft_dynamics_slow!

function spacecraft_dynamics_fast_control!(du::ComponentVector, u::ComponentVector, p, t::Float64)
    sc_state = u.sc
    sc_du = du.sc
    spacecraft = p.args.dynamics_model.spacecraft
    debug_control = p.shared_buffers.debug_control[]
    p.shared_buffers.current_time[] = t
    use_rhs_batch = _rhs_batch_parallel_enabled(length(spacecraft))
    if use_rhs_batch
        minbatch = max(1, Int(ceil(length(spacecraft) / Polyester.num_cores())))
        @batch minbatch=minbatch for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                du_view.pos .= 0.0
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.0
                    du_view.ω .= inertia_tensor \ τ_body
                end

                du_view.heat_loads .= 0.0
            end
        end
    else
        @inbounds for i in eachindex(sc_state)
            if !p.is_active[i]
                sc_du[i] .= 0.0
                continue
            end
            @views begin
                sc_view = sc_state[i]
                du_view = sc_du[i]
                forces = MVector{3, Float64}(0.0, 0.0, 0.0)
                torques = MVector{3, Float64}(0.0, 0.0, 0.0)
                mass_rate = _accumulate_control_effectors!(forces, torques, sc_view, p, i, t, debug_control)

                du_view.pos .= 0.0
                du_view.vel .= forces / sc_view.mass
                du_view.mass = mass_rate

                if p.args.mission_configuration.orientation_sim
                    inertia_tensor = spacecraft[i].inertia_tensor
                    τ_body = SVector{3, Float64}(torques)
                    du_view.q .= 0.0
                    du_view.ω .= inertia_tensor \ τ_body
                end

                du_view.heat_loads .= 0.0
            end
        end
    end
end # function spacecraft_dynamics_fast_control!

function build_initial_conditions(args)::ComponentVector
    # 1. Build the structure (Axis) based on each spacecraft's unique body count
    # This identifies exactly how many heat_load slots each SC needs
    sc_shapes = map(args.dynamics_model.spacecraft) do sc
        # Get the number of bodies for this specific spacecraft
        n_bodies = length(sc.links)
        mass = sc.dry_mass + sc.prop_mass
        # Create the initial state for this spacecraft with the correct size for heat_loads
        if args.mission_configuration.orientation_sim
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies), # Variable size!
                q = Float64[0.0, 0.0, 0.0, 1.0], 
                ω = zeros(3)
            )
        else
            return (
                pos = zeros(3), 
                vel = zeros(3), 
                mass = mass, 
                heat_loads = zeros(n_bodies)
            )
        end
    end

    # 2. Pack everything into one ComponentVector
    # Julia allocates ONE flat array and calculates all offsets automatically
    state = ComponentVector(sc = sc_shapes) # Add more components here as needed in the future (e.g., separate orientation state if not using quaternions, etc.)

    # 3. Fill the values (The logic remains the same)
    for i in eachindex(args.dynamics_model.spacecraft)
        spacecraft = args.dynamics_model.spacecraft[i]
        sc_view = state.sc[i]
        r0, v0 = orbitalelemtorv(spacecraft.initial_condition, args.environment_model.planet)
        sc_view.pos .= r0
        sc_view.vel .= v0
        # sc_view.mass .= spacecraft.dry_mass + spacecraft.prop_mass
        # Note: heat_loads is already the correct size for this specific i!
        sc_view.heat_loads .= 0.0  
        
        if args.mission_configuration.orientation_sim
            sc_view.q .= spacecraft.initial_condition.q
            sc_view.ω .= spacecraft.initial_condition.ang_vel
        end
    end

    return state
end
