include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))
using LinearAlgebra
using StaticArrays
using LoopVectorization
using ComponentArrays
using OrdinaryDiffEq
using DiffEqCallbacks
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
const NBODY_EPHEMERIS_CACHE_SCHEMA_VERSION = "1"
const _EPHEMERIS_REUSE_LOCK = ReentrantLock()
const SRPEphemerisReuseKey = Tuple{String, Int64, Int64, Int64}
const NBodyEphemerisReuseKey = Tuple{String, String, Int64, Int64, Int64}
const PlanetFrameEphemerisReuseKey = Tuple{String, String, NTuple{9, Int64}, Int64, Int64, Int64}
const _SRP_EPHEMERIS_REUSE_CACHE = Dict{SRPEphemerisReuseKey, SimulationModel.SRPSunEphemerisCache}()
const _NBODY_EPHEMERIS_REUSE_CACHE = Dict{NBodyEphemerisReuseKey, SimulationModel.NBodyEphemerisCache}()
const _PLANET_FRAME_EPHEMERIS_REUSE_CACHE = Dict{PlanetFrameEphemerisReuseKey, SimulationModel.PlanetFrameEphemerisCache}()
const _NBODY_EPHEMERIS_PREWARMED_CACHE = Dict{NBodyEphemerisReuseKey, SimulationModel.NBodyEphemerisCache}()

@inline _typed_normalize_warning_enabled() = _engine_env_get("SPACEAGORA_WARN_NORMALIZE", "1") == "1"
@inline _typed_allow_transition_normalize() = _engine_env_get("SPACEAGORA_ALLOW_TYPED_NORMALIZE", "0") == "1"
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
    return (
        round(Int64, planet.ω[1] * 1e12),
        round(Int64, planet.ω[2] * 1e12),
        round(Int64, planet.ω[3] * 1e12),
        round(Int64, planet.α * 1e12),
        round(Int64, planet.δ * 1e12),
        round(Int64, planet.Rp_e * 1e6),
        round(Int64, planet.Rp_p * 1e6),
        round(Int64, planet.Rp_m * 1e6),
        round(Int64, planet.μ * 1e-6),
    )
end

@inline function _body_query_names_reuse_key(body_query_names::AbstractVector{String})::String
    return join(body_query_names, '\0')
end

@inline function _ephemerides_model_reuse_key(ephemerides_model)::String
    cache_key = SimulationModel.ephemerides_cache_key(ephemerides_model)
    if cache_key == (:spice,)
        return "spice"
    end
    return "simple:$(cache_key[2]):$(cache_key[3])"
end

@inline function _srp_ephemeris_reuse_key(primary_body_name::String, et_start::Float64, mission_end_s::Float64, dt_s::Float64)::SRPEphemerisReuseKey
    return (
        primary_body_name,
        _cache_time_key(et_start),
        _cache_time_key(mission_end_s),
        _cache_time_key(dt_s)
    )
end

@inline function _nbody_ephemeris_reuse_key(primary_body_name::String, body_query_names::Vector{String}, et_start::Float64, mission_end_s::Float64, dt_s::Float64)::NBodyEphemerisReuseKey
    return (
        primary_body_name,
        _body_query_names_reuse_key(body_query_names),
        _cache_time_key(et_start),
        _cache_time_key(mission_end_s),
        _cache_time_key(dt_s)
    )
end

@inline function _planet_frame_ephemeris_reuse_key(planet, ephemerides_model, et_start::Float64, mission_end_s::Float64, dt_s::Float64)::PlanetFrameEphemerisReuseKey
    return (
        string(planet.name),
        _ephemerides_model_reuse_key(ephemerides_model),
        _planet_transform_key(planet),
        _cache_time_key(et_start),
        _cache_time_key(mission_end_s),
        _cache_time_key(dt_s)
    )
end

@inline function _ephemeris_reuse_lookup(cache::AbstractDict{K, T}, key::K) where {K, T}
    return lock(_EPHEMERIS_REUSE_LOCK) do
        get(cache, key, nothing)
    end
end

@inline function _ephemeris_reuse_store!(cache::AbstractDict{K, T}, key::K, value::T, max_entries::Int)::T where {K, T}
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

@inline function _ephemeris_explicit_cache_store!(cache::AbstractDict{K, T}, key::K, value::T; replace::Bool=false)::T where {K, T}
    return lock(_EPHEMERIS_REUSE_LOCK) do
        if !replace
            existing = get(cache, key, nothing)
            if existing !== nothing
                return existing
            end
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
        empty!(_NBODY_EPHEMERIS_PREWARMED_CACHE)
    end
    return nothing
end

@inline function _has_active_srp_effector(dynamic_effectors::Tuple)::Bool
    @inbounds for effector in dynamic_effectors
        if effector isa SimulationModel.SolarRadiationPressureModel &&
           effector.A > 0.0 &&
           (effector.direct || effector.albedo)
            return true
        end
    end
    return false
end

@inline function _is_nbody_effector_like(effector)::Bool
    return effector isa SimulationModel.NBodyGravityModel ||
           (hasproperty(effector, :body_names) && hasproperty(effector, :primary_body_name))
end

@inline function _effector_parallel_mode()::Symbol
    return SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_EFFECTOR_PARALLEL"; default="auto")
end

@inline function _rhs_batch_parallel_mode()::Symbol
    return SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_RHS_BATCH_PARALLEL"; default="auto")
end

@inline function _rhs_batch_thread_threshold()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_BATCH_THREAD_THRESHOLD", 16)
end

@inline function _profile_forces_serial_rhs()::Bool
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_PARALLEL_PROFILE", "")))
    return raw in ("r0", "serial", "r0_true_serial", "true_serial")
end

@inline function _rhs_batch_parallel_enabled(env::SimulationModel.RhsPlanEnvConfig, num_spacecraft::Int)::Bool
    if env.profile_forces_serial
        return false
    end
    mode = env.batch_parallel_mode
    if mode == :off
        return false
    elseif mode == :on
        return true
    end
    return num_spacecraft >= env.batch_thread_threshold && Polyester.num_cores() > 1
end

@inline function _rhs_batch_parallel_enabled(p, num_spacecraft::Int)::Bool
    return _rhs_batch_parallel_enabled(_rhs_env_config(p), num_spacecraft)
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

@inline function _rhs_flat_packet_target_min_ns()::Float64
    return _parse_positive_float_env("SPACEAGORA_RHS_FLAT_PACKET_TARGET_MIN_NS", 2.5e4)
end

@inline function _rhs_flat_packet_scheduler_mode()::Symbol
    return SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_RHS_FLAT_PACKET_SCHEDULER"; default="auto")
end

@inline function _rhs_flat_packet_min_items()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_FLAT_PACKET_MIN_ITEMS", 128)
end

@inline function _rhs_flat_packet_work_ns_threshold()::Float64
    return _parse_positive_float_env("SPACEAGORA_RHS_FLAT_PACKET_WORK_NS_THRESHOLD", 5.0e6)
end

@inline function _rhs_flat_packet_heterogeneity_threshold()::Float64
    return _parse_positive_float_env("SPACEAGORA_RHS_FLAT_PACKET_HETEROGENEITY_THRESHOLD", 3.0)
end

@inline function _rhs_flat_packet_overhead_disable_ratio()::Float64
    return _parse_unit_float_env("SPACEAGORA_RHS_FLAT_PACKET_OVERHEAD_DISABLE_RATIO", 0.10)
end

@inline function _rhs_flat_packet_overhead_min_samples()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_FLAT_PACKET_OVERHEAD_MIN_SAMPLES", 4)
end

@inline function _rhs_effector_cost_min_samples()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_EFFECTOR_COST_MIN_SAMPLES", 4)
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

@inline function _effector_observed_cost_ns_per_item(env::SimulationModel.RhsPlanEnvConfig, shared_buffers)::Float64
    default_cost = env.effector_cost_ns_per_item_default
    shared_buffers === nothing && return default_cost
    if !(hasproperty(shared_buffers, :effector_cost_ns_per_item) && hasproperty(shared_buffers, :effector_cost_samples))
        return default_cost
    end
    samples = Int(getproperty(shared_buffers, :effector_cost_samples)[])
    estimate = Float64(getproperty(shared_buffers, :effector_cost_ns_per_item)[])
    if samples >= env.effector_cost_min_samples && isfinite(estimate) && estimate > 0.0
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
    α = _rhs_env_config_from_buffers(shared_buffers).effector_cost_ema_alpha
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

function _ensure_rhs_effector_cost_model!(shared_buffers, n_effectors::Int)::Nothing
    shared_buffers === nothing && return nothing
    n_effectors >= 0 || return nothing
    if hasproperty(shared_buffers, :rhs_effector_cost_ns)
        costs = shared_buffers.rhs_effector_cost_ns[]
        if length(costs) < n_effectors
            old_len = length(costs)
            resize!(costs, n_effectors)
            @inbounds for idx in (old_len + 1):n_effectors
                costs[idx] = NaN
            end
        end
    end
    if hasproperty(shared_buffers, :rhs_effector_cost_samples)
        samples = shared_buffers.rhs_effector_cost_samples[]
        if length(samples) < n_effectors
            old_len = length(samples)
            resize!(samples, n_effectors)
            @inbounds for idx in (old_len + 1):n_effectors
                samples[idx] = Int64(0)
            end
        end
    end
    return nothing
end

@inline function _rhs_effector_observed_cost_ns(shared_buffers, eff_idx::Int, fallback_ns::Float64)::Float64
    shared_buffers === nothing && return fallback_ns
    if !(hasproperty(shared_buffers, :rhs_effector_cost_ns) && hasproperty(shared_buffers, :rhs_effector_cost_samples))
        return fallback_ns
    end
    costs = shared_buffers.rhs_effector_cost_ns[]
    samples = shared_buffers.rhs_effector_cost_samples[]
    if eff_idx <= length(costs) && eff_idx <= length(samples)
        estimate = Float64(costs[eff_idx])
        if samples[eff_idx] >= _rhs_env_config_from_buffers(shared_buffers).rhs_effector_cost_min_samples && isfinite(estimate) && estimate > 0.0
            return estimate
        end
    end
    return fallback_ns
end

function _update_rhs_effector_cost_model!(
    shared_buffers,
    eff_idx::Int,
    elapsed_ns::Float64,
)::Nothing
    shared_buffers === nothing && return nothing
    elapsed_ns > 0.0 || return nothing
    if !(hasproperty(shared_buffers, :rhs_effector_cost_ns) && hasproperty(shared_buffers, :rhs_effector_cost_samples))
        return nothing
    end
    _ensure_rhs_effector_cost_model!(shared_buffers, eff_idx)
    costs = shared_buffers.rhs_effector_cost_ns[]
    samples = shared_buffers.rhs_effector_cost_samples[]
    α = _rhs_env_config_from_buffers(shared_buffers).effector_cost_ema_alpha
    previous = costs[eff_idx]
    costs[eff_idx] = if isfinite(previous) && previous > 0.0
        (1.0 - α) * previous + α * elapsed_ns
    else
        elapsed_ns
    end
    samples[eff_idx] = min(typemax(Int64), samples[eff_idx] + Int64(1))
    return nothing
end

@inline function _dynamic_effector_thread_decision(
    args::SimulationConfiguration,
    p,
    dynamic_effectors::Tuple,
    num_sats::Int
)
    return _dynamic_effector_thread_decision(
        _rhs_env_config(p), _policy_env_config(p), args, p, dynamic_effectors, num_sats
    )
end

@inline function _dynamic_effector_thread_decision(
    env::SimulationModel.RhsPlanEnvConfig,
    penv::Union{Nothing, SimulationModel.PolicyDecisionEnvConfig},
    args::SimulationConfiguration,
    p,
    dynamic_effectors::Tuple,
    num_sats::Int
)
    mode = env.effector_parallel_mode
    n_effectors = length(dynamic_effectors)
    if n_effectors <= 1
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    if !_dynamic_effectors_parallel_supported(dynamic_effectors)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end

    outer_active = penv === nothing ? _effector_outer_parallel_hint() : penv.outer_parallel_active
    allow_with_outer = env.effector_allow_with_outer
    budget = penv === nothing ?
        SimulationModel.ParallelPolicy.effective_inner_thread_budget() : penv.inner_thread_budget
    share_budget = _effector_satellite_share_budget(num_sats, budget)
    if outer_active && allow_with_outer
        share_budget = max(1, fld(share_budget, 2))
    end
    inner_floor = (!outer_active && num_sats > 1 && budget > 1) ? min(2, budget) : 1
    max_allotment = min(env.effector_max_threads, budget, max(share_budget, inner_floor))

    shared_buffers = _effector_shared_buffers(p)
    per_effector_cost_ns = _effector_observed_cost_ns_per_item(env, shared_buffers)
    estimated_work_ns = per_effector_cost_ns * n_effectors
    work_per_worker_ns = estimated_work_ns / max(1, max_allotment)
    target_ns = env.effector_work_ns_per_worker_threshold * (outer_active ? env.effector_outer_work_scale : 1.0)
    heavy_work = work_per_worker_ns >= target_ns

    policy = SimulationModel.ParallelPolicy.thread_policy_decision(
        n_effectors;
        mode=mode,
        threshold=env.effector_thread_threshold,
        heavy_work=heavy_work,
        heavy_only=env.effector_heavy_only,
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:dynamic_effectors,
        env=penv
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

# Returns true when satellite-level batching occupies every available thread,
# making nested per-effector threading purely additive overhead.
@inline function _satellite_batch_saturates_pool(active_sats::Int, budget::Int)::Bool
    return active_sats >= budget && budget > 1
end

# Effectors whose flat-mode routing benefit satellite_batch cannot replicate at any
# thread count: harmonics gets the cache-friendly SIMD batch kernel across the whole
# satellite batch (satellite_batch's Polyester loop instead calls the ordinary
# per-satellite Pines recursion once per satellite), and (J2-)inverse-square gravity
# gets a batched pre-pass specifically to avoid per-satellite Polyester dispatch
# overhead for a kernel too cheap to pay for it (_accumulate_invsq_flat_batch! /
# _accumulate_invsq_j2_flat_batch! in dynamics_rhs.jl). NBody/SRP are deliberately not
# included: satellite_batch already hoists their expensive shared ephemeris lookups
# via _prefill_shared_body_samples!, so it doesn't lose anything structural for those.
#
# This matters because n_effectors >= env.flat_min_effectors (default 3) below is
# gating flat-queue eligibility on "enough effector *types* to justify the queue",
# which has nothing to do with whether one of those effectors has an algorithmic
# benefit independent of effector count. A 2-effector harmonics+drag constellation
# (a common shape) would otherwise never reach the flat queue and would fall through
# to _satellite_batch_saturates_pool once active_sats >= budget, silently discarding
# the harmonics batch kernel's win even though it has nothing to do with effector count.
@inline function _rhs_flat_batch_privileged_effector(effector)::Bool
    effector isa SimulationModel.GravitationalHarmonicsModel && return true
    effector isa SimulationModel.InverseSquaredGravityModel && return !effector.gravity_gradient
    effector isa SimulationModel.InverseSquaredJ2GravityModel && return !effector.gravity_gradient
    return false
end

@inline function _rhs_flat_has_batch_privileged_effector(dynamic_effectors::Tuple)::Bool
    @inbounds for effector in dynamic_effectors
        _rhs_flat_batch_privileged_effector(effector) && return true
    end
    return false
end

# Forces an effector decision to serial, preserving the mode/policy fields for
# telemetry while making it structurally impossible to enable nested effector
# threads under satellite_batch.
@inline function _with_serial_effector_decision(effector_decision)
    return (
        use_threads=false,
        allotment=1,
        mode=effector_decision.mode,
        policy_applied=effector_decision.policy_applied,
    )
end

@inline function _rhs_execution_mode_env()::Symbol
    raw = lowercase(strip(_engine_env_get("SPACEAGORA_RHS_EXECUTION_MODE", "auto")))
    if raw in ("auto", "")
        return :auto
    elseif raw in ("serial", "off", "none")
        return :serial
    elseif raw in ("satellite", "satellite_batch", "batch")
        return :satellite_batch
    elseif raw in ("per_satellite", "per_satellite_effector_reduce", "effector")
        return :per_satellite_effector_reduce
    elseif raw in ("flat", "flat_constellation", "flat_constellation_effector_queue")
        return :flat_constellation_effector_queue
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_RHS_EXECUTION_MODE='$raw'. Use one of: auto, serial, satellite, per_satellite, flat."
    ))
end

@inline function _rhs_flat_min_sats()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_FLAT_MIN_SATS", 24)
end

@inline function _rhs_flat_min_effectors()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_FLAT_MIN_EFFECTORS", 3)
end

@inline function _rhs_flat_work_ns_threshold()::Float64
    return _parse_positive_float_env("SPACEAGORA_EFFECTOR_FLAT_WORK_NS_THRESHOLD", 2.0e5)
end

@inline function _rhs_flat_work_per_worker_ns_threshold()::Float64
    return _parse_positive_float_env("SPACEAGORA_RHS_FLAT_WORK_PER_WORKER_NS_THRESHOLD", 2.5e4)
end

@inline function _rhs_flat_cost_heterogeneity_threshold()::Float64
    return _parse_positive_float_env("SPACEAGORA_EFFECTOR_FLAT_COST_HETEROGENEITY_THRESHOLD", 3.0)
end

@inline function _rhs_flat_min_thread_budget()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env(
        "SPACEAGORA_EFFECTOR_FLAT_MIN_THREAD_BUDGET",
        SimulationModel.ParallelPolicy.auto_thread_min_budget(),
    )
end

@inline function _rhs_harmonics_flat_experimental_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_HARMONICS_FLAT_EXPERIMENTAL", false)
end

@inline function _rhs_harmonics_batch_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_HARMONICS_BATCH_ENABLED", true) ||
        _rhs_harmonics_flat_experimental_enabled()
end

# Minimum satellites per worker for the harmonics SIMD batch (SIMD efficiency floor).
# The formula workers = min(budget, fld(active_sats, floor)) automatically scales with
# both satellite count and thread budget.  Default 4: one AVX2 / half AVX-512 register
# width — enough for @turbo to have a full SIMD iteration while allowing 16 workers at
# 64 satellites (fld(64,4)=16) and 32 workers at 128 satellites (fld(128,4)=32).
# Raise to 8 if targeting AVX-512-only machines where 4-wide @turbo is suboptimal.
@inline function _rhs_harmonics_batch_min_sats_per_worker()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env(
        "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER", 4
    )
end

# Unlike _dynamic_effector_thread_decision and the density/control/thermal
# callback paths, the harmonics-batch flat-constellation route did not
# previously consult outer_parallel_active() at all -- it fires on every
# RHS/ODE step (not once per sample), so an outer worker already blocked in
# Threads.@sync/Base.@sync repeatedly spawned its own nested batch of workers
# throughout the whole integration. That's the mechanism behind the severe,
# livelock-like nested outer+inner contention documented in
# THREAD_ALLOCATION_AND_GRAM_CONCURRENCY_HANDOFF.md Finding 3 (some points
# ranged from ~1x overhead to a multi-minute hang on identical repeated
# runs). Default false, matching _effector_allow_with_outer's default and
# this codebase's own documented recommendation to never split the thread
# budget between outer and inner parallelism.
@inline function _harmonics_batch_allow_with_outer()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER", false)
end

"""
    _snapshot_rhs_plan_env_config() -> SimulationModel.RhsPlanEnvConfig

Resolve every env-derived knob consulted by the per-RHS-call execution-plan
routing chain into a typed snapshot.  Built once at run_simulation setup
(inside any active SimulationEngineConfig override scope); hot paths read
plain struct fields via `_rhs_env_config(p)` instead of re-parsing ENV per
RHS evaluation.
"""
function _snapshot_rhs_plan_env_config()::SimulationModel.RhsPlanEnvConfig
    return SimulationModel.RhsPlanEnvConfig(
        _rhs_execution_mode_env(),
        _profile_forces_serial_rhs(),
        _rhs_batch_parallel_mode(),
        _rhs_batch_thread_threshold(),
        _effector_parallel_mode(),
        _effector_thread_threshold(),
        _effector_max_threads(),
        _effector_allow_with_outer(),
        _effector_heavy_only(),
        _effector_cost_ns_per_item_default(),
        _effector_cost_min_samples(),
        _effector_cost_ema_alpha(),
        _effector_work_ns_per_worker_threshold(),
        _effector_outer_work_scale(),
        _rhs_flat_min_sats(),
        _rhs_flat_min_effectors(),
        _rhs_flat_work_ns_threshold(),
        _rhs_flat_work_per_worker_ns_threshold(),
        _rhs_flat_cost_heterogeneity_threshold(),
        _rhs_flat_min_thread_budget(),
        _rhs_harmonics_batch_enabled(),
        _rhs_harmonics_batch_min_sats_per_worker(),
        SimulationModel.ParallelPolicy.harmonics_batch_spin_barrier_enabled(),
        _harmonics_batch_allow_with_outer(),
        _rhs_effector_cost_min_samples(),
        _rhs_flat_packet_target_min_ns(),
        _rhs_flat_packet_scheduler_mode(),
        _rhs_flat_packet_min_items(),
        _rhs_flat_packet_work_ns_threshold(),
        _rhs_flat_packet_heterogeneity_threshold(),
        _rhs_flat_packet_overhead_disable_ratio(),
        _rhs_flat_packet_overhead_min_samples(),
    )
end

# Run-scoped snapshot accessors.  Fall back to live ENV parsing when the
# snapshot is unset (hand-constructed ODEParams in unit tests / withenv
# probes) so read-ENV-at-use behavior is preserved outside run_simulation.
@inline function _rhs_env_config_from_buffers(shared_buffers)::SimulationModel.RhsPlanEnvConfig
    if shared_buffers !== nothing && hasproperty(shared_buffers, :rhs_env_config)
        cfg = shared_buffers.rhs_env_config[]
        cfg === nothing || return cfg
    end
    return _snapshot_rhs_plan_env_config()
end

@inline _rhs_env_config(p)::SimulationModel.RhsPlanEnvConfig =
    _rhs_env_config_from_buffers(_effector_shared_buffers(p))

@inline function _policy_env_config(p)::Union{Nothing, SimulationModel.PolicyDecisionEnvConfig}
    sb = _effector_shared_buffers(p)
    if sb !== nothing && hasproperty(sb, :policy_env_config)
        return sb.policy_env_config[]
    end
    return nothing
end

# Resolve all run-scoped env snapshots (policy, RHS plan, callbacks) onto the
# shared buffers.  Called once from run_simulation setup so every value is
# captured inside the active SimulationEngineConfig override scope; the RHS
# and callback hot paths then read plain struct fields.
function _initialize_runtime_env_config!(p)
    p.shared_buffers.policy_env_config[] = SimulationModel.ParallelPolicy.snapshot_policy_decision_env()
    p.shared_buffers.rhs_env_config[] = _snapshot_rhs_plan_env_config()
    p.shared_buffers.callback_env_config[] = SimulationModel.SimulationCallbacks._snapshot_callback_env_config()
    return nothing
end

@inline function _rhs_flat_supported(dynamic_effectors::Tuple)::Bool
    if length(dynamic_effectors) == 1
        return (dynamic_effectors[1] isa SimulationModel.GravitationalHarmonicsModel &&
                _dynamic_effector_threadsafe(dynamic_effectors[1]) &&
                _rhs_harmonics_batch_enabled()) ||
            _rhs_single_invsq_flat_supported(dynamic_effectors)
    end
    return length(dynamic_effectors) > 1 && _dynamic_effectors_parallel_supported(dynamic_effectors)
end

@inline function _rhs_flat_supported(env::SimulationModel.RhsPlanEnvConfig, dynamic_effectors::Tuple)::Bool
    if length(dynamic_effectors) == 1
        return (dynamic_effectors[1] isa SimulationModel.GravitationalHarmonicsModel &&
                _dynamic_effector_threadsafe(dynamic_effectors[1]) &&
                env.harmonics_batch_enabled) ||
            _rhs_single_invsq_flat_supported(dynamic_effectors)
    end
    return length(dynamic_effectors) > 1 && _dynamic_effectors_parallel_supported(dynamic_effectors)
end

@inline function _rhs_single_harmonics_flat_supported(env::SimulationModel.RhsPlanEnvConfig, dynamic_effectors::Tuple)::Bool
    return length(dynamic_effectors) == 1 &&
        dynamic_effectors[1] isa SimulationModel.GravitationalHarmonicsModel &&
        _dynamic_effector_threadsafe(dynamic_effectors[1]) &&
        env.harmonics_batch_enabled
end

# Single-effector fast path for plain (J2-)inverse-square gravity, mirroring
# _rhs_single_harmonics_flat_supported. Unlike harmonics, there is no coefficient
# table to batch over — the payoff is routing around per-satellite Polyester
# dispatch entirely in favor of the serial batchable-effector pre-pass
# (_accumulate_invsq_flat_batch!/_accumulate_invsq_j2_flat_batch!), whose per-item
# cost is too small to ever amortise task-spawn overhead. Excluded when
# `gravity_gradient=true` since the batch kernel writes force only.
@inline function _rhs_single_invsq_flat_supported(dynamic_effectors::Tuple)::Bool
    length(dynamic_effectors) == 1 || return false
    effector = dynamic_effectors[1]
    return (effector isa SimulationModel.InverseSquaredGravityModel ||
            effector isa SimulationModel.InverseSquaredJ2GravityModel) &&
        !effector.gravity_gradient &&
        _dynamic_effector_threadsafe(effector)
end

@inline function _rhs_invsq_flat_min_sats()::Int
    return SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_INVSQ_FLAT_MIN_SATS", 8)
end

@inline function _rhs_effectors_have_heavy_or_heterogeneous_cost(
    dynamic_effectors::Tuple,
    heterogeneity_threshold::Float64
)::Bool
    has_nbody = false
    has_aero = false
    has_harmonics = false
    min_cost = Inf
    max_cost = 0.0
    @inbounds for effector in dynamic_effectors
        cost = _rhs_effector_static_cost_ns(effector)
        min_cost = min(min_cost, cost)
        max_cost = max(max_cost, cost)
        if effector isa SimulationModel.NBodyGravityModel
            has_nbody = true
        elseif effector isa SimulationModel.AerodynamicCoefficientfM
            has_aero = true
        elseif effector isa SimulationModel.GravitationalHarmonicsModel
            has_harmonics = true
        end
    end
    heterogeneity = max_cost / max(min_cost, 1.0)
    return has_nbody || has_aero || has_harmonics ||
        heterogeneity >= heterogeneity_threshold
end

# Default off: this is a new, not-yet-broadly-validated optimization (see
# rhs_plan_step_cache's docstring in runtime_types.jl for the mechanism). Every
# other perf-sensitive toggle in this file defaults to today's proven behavior
# and requires an explicit opt-in; this follows the same convention.
@inline function _rhs_plan_step_cache_enabled()::Bool
    return SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_RHS_PLAN_STEP_CACHE", false)
end

# Tried and reverted: routing _spacecraft_dynamics_flat_constellation_effector_queue!'s
# final per-satellite state-assembly pass (dynamics_rhs.jl) through
# threaded_foreach_persistent instead of plain threaded_foreach (which spawns a fresh
# Threads.@spawn/Threads.@sync task batch on every RHS call). Motivated by profiling a
# no_gram/512-sat/monolithic scenario: Julia's own multiq/work-stealing scheduler
# internals (multiq_deletemin, julia_multiq_check_empty) become a top-20 self-time
# cost at 8 threads that's absent at 4, and thread 1's utilization jumps from 28% to
# 64% between 4 and 8 threads -- well before any GRAM/lock/density-callback code is in
# the picture. Measured directly (SPACEAGORA_RHS_FLAT_ASSEMBLY_PERSISTENT_DISPATCH,
# same thread ladder, same scenario): correctness was clean (full test suite +
# trajectory parity), but the persistent-pool path was slower, not faster --
# 0.9x at 2 threads, 0.811x at 4 threads (worse at the previously-best point), ~1.0x
# at 8 (no help at all). This assembly loop's per-item work (assigning already-
# computed forces/torques into du) is cheap enough that the persistent pool's own
# dispatch overhead (locked pool lookup, channel-based wake signaling) costs more
# than plain Threads.@spawn for a payload this light. Root cause of the original
# 4-vs-8-thread regression is still open; this specific fix doesn't work.

# Thin caching wrapper around _rhs_execution_plan_uncached: when enabled, reuses
# the same routing decision across every solver stage within one accepted step
# instead of re-deriving it (active-satellite count, effector cost decision,
# routing heuristics) on every single stage call. See rhs_plan_step_cache's
# docstring (runtime_types.jl) for how/when the cache is invalidated.
@inline function _rhs_execution_plan(
    args::SimulationConfiguration,
    p,
    dynamic_effectors::Tuple,
    num_sats::Int
)::SimulationModel.RhsExecutionPlan
    if p !== nothing && hasproperty(p, :shared_buffers) && _rhs_plan_step_cache_enabled()
        cached = p.shared_buffers.rhs_plan_step_cache[]
        cached === nothing || return cached
        plan = _rhs_execution_plan_uncached(args, p, dynamic_effectors, num_sats)
        p.shared_buffers.rhs_plan_step_cache[] = plan
        return plan
    end
    return _rhs_execution_plan_uncached(args, p, dynamic_effectors, num_sats)
end

@inline function _rhs_execution_plan_uncached(
    args::SimulationConfiguration,
    p,
    dynamic_effectors::Tuple,
    num_sats::Int
)::SimulationModel.RhsExecutionPlan
    # Calibration override: return the pre-measured best plan without heuristic routing.
    if p !== nothing && hasproperty(p, :shared_buffers)
        override = p.shared_buffers.rhs_plan_override[]
        override === nothing || return override
    end
    env = _rhs_env_config(p)
    penv = _policy_env_config(p)
    forced_mode = env.execution_mode
    n_effectors = length(dynamic_effectors)
    active_sats = if p === nothing || !hasproperty(p, :is_active)
        num_sats
    else
        count(identity, p.is_active)
    end
    budget = penv === nothing ?
        SimulationModel.ParallelPolicy.effective_inner_thread_budget() : penv.inner_thread_budget
    outer_active = penv === nothing ?
        SimulationModel.ParallelPolicy.outer_parallel_active() : penv.outer_parallel_active
    # Compute the threading preference from the policy; route selection below may
    # override it via _with_serial_effector_decision when satellite_batch is chosen.
    effector_decision = _dynamic_effector_thread_decision(env, penv, args, p, dynamic_effectors, num_sats)

    # ── Forced modes ────────────────────────────────────────────────────────────
    # Satellite_batch and serial always disable inner effector threading:
    # the outer @batch loop already saturates available threads.
    if forced_mode == :serial
        return (
            mode=:serial,
            allotment=1,
            scheduler=:serial,
            dominant_axis=:serial,
            policy_applied=true,
            effector_decision=_with_serial_effector_decision(effector_decision),
        )
    elseif forced_mode == :satellite_batch
        return (
            mode=:satellite_batch,
            allotment=1,
            scheduler=:static,
            dominant_axis=:satellite,
            policy_applied=true,
            effector_decision=_with_serial_effector_decision(effector_decision),
        )
    elseif forced_mode == :per_satellite_effector_reduce
        return (
            mode=:per_satellite_effector_reduce,
            allotment=effector_decision.allotment,
            scheduler=:auto,
            dominant_axis=:per_satellite_inner_effector,
            policy_applied=true,
            effector_decision=effector_decision,
        )
    elseif forced_mode == :flat_constellation_effector_queue
        # Same nested-outer-split hazard as the auto-routed flat/harmonics
        # branches below -- an explicit SPACEAGORA_RHS_EXECUTION_MODE=flat
        # request is still subject to it, since this route fires every RHS
        # step and the other outer_parallel_active-aware call sites in this
        # codebase (density/control/thermal callbacks, dynamic effectors)
        # don't distinguish forced vs. auto requests either.
        if active_sats <= 1 || !_rhs_flat_supported(env, dynamic_effectors) || budget <= 1 ||
           (outer_active && !env.harmonics_batch_allow_with_outer)
            return (
                mode=:satellite_batch,
                allotment=1,
                scheduler=:static,
                dominant_axis=:satellite,
                policy_applied=true,
                effector_decision=_with_serial_effector_decision(effector_decision),
            )
        end
        return (
            mode=:flat_constellation_effector_queue,
            allotment=min(max(1, budget), max(active_sats * n_effectors, 1)),
            scheduler=:dynamic,
            dominant_axis=:flat_effector,
            policy_applied=true,
            effector_decision=_with_serial_effector_decision(effector_decision),
        )
    end

    # ── Auto routing ─────────────────────────────────────────────────────────────

    # Single-harmonics SIMD batch: route through flat_constellation whenever the
    # constellation is large enough to give each worker at least min_sats_per_worker
    # satellites (viable_workers >= 2), or when running single-threaded (budget <= 1).
    # The SIMD batch kernel loads each spherical harmonic coefficient pair once and
    # broadcasts it to the full satellite batch via @turbo, which is more cache-efficient
    # than Polyester @batch where every thread independently traverses the coefficient
    # table for its satellite slice.  The worker count is capped at viable_workers inside
    # _accumulate_harmonics_flat_batch!, so this routing is safe at any thread budget.
    single_harmonics_flat = _rhs_single_harmonics_flat_supported(env, dynamic_effectors)
    if single_harmonics_flat && active_sats >= env.flat_min_sats && active_sats > 1
        # A higher-level campaign (or benchmark harness) already owns an outer
        # split -- this route fires on every RHS/ODE step, not once per
        # sample, so nesting its own multi-worker batch here would repeatedly
        # oversubscribe the same thread pool an already-blocked outer worker
        # is waiting on, throughout the whole integration (see
        # _harmonics_batch_allow_with_outer's docstring / Finding 3). Force
        # serial instead, same as this function's other serial fallbacks.
        if outer_active && !env.harmonics_batch_allow_with_outer
            return (
                mode=:satellite_batch,
                allotment=1,
                scheduler=:static,
                dominant_axis=:satellite,
                policy_applied=true,
                effector_decision=_with_serial_effector_decision(effector_decision),
            )
        end
        min_sats_floor = env.harmonics_batch_spin_barrier ?
            1 : env.harmonics_batch_min_sats_per_worker
        viable_workers = fld(active_sats, max(1, min_sats_floor))
        if budget <= 1 || viable_workers >= 2
            return (
                mode=:flat_constellation_effector_queue,
                allotment=min(max(1, budget), viable_workers),
                scheduler=:dynamic,
                dominant_axis=:flat_effector,
                policy_applied=true,
                effector_decision=_with_serial_effector_decision(effector_decision),
            )
        end
    end

    # Single inverse-square (J2-)gravity fast path: the effector body is a few
    # FLOPs, far too cheap to ever amortise Polyester per-satellite task-spawn
    # overhead (unlike harmonics, there's no per-worker SIMD batch to size —
    # the win is entirely from replacing satellite_batch's per-task effector
    # dispatch with a single serial pre-pass). Threading is reserved for the
    # subsequent per-satellite RHS-assembly pass, which every route (including
    # satellite_batch) already parallelizes.
    single_invsq_flat = _rhs_single_invsq_flat_supported(dynamic_effectors)
    if single_invsq_flat && active_sats >= _rhs_invsq_flat_min_sats() && active_sats > 1
        # Same nested-outer-split hazard as the harmonics-batch route above --
        # this also fires every RHS step and shares the same multi-worker flat
        # batch kernel.
        if outer_active && !env.harmonics_batch_allow_with_outer
            return (
                mode=:satellite_batch,
                allotment=1,
                scheduler=:static,
                dominant_axis=:satellite,
                policy_applied=true,
                effector_decision=_with_serial_effector_decision(effector_decision),
            )
        end
        return (
            mode=:flat_constellation_effector_queue,
            allotment=min(max(1, budget), active_sats),
            scheduler=:dynamic,
            dominant_axis=:flat_effector,
            policy_applied=true,
            effector_decision=_with_serial_effector_decision(effector_decision),
        )
    end

    # Safety: not enough satellites, threads, or thread-safe effectors for any
    # satellite-parallel path. For single-satellite cases with a thread budget,
    # the full budget is available for inner effector parallelism — preserve
    # the effector_decision rather than forcing it serial. For the budget<=1 or
    # non-flat-supported cases, serial is correct on both axes.
    if active_sats <= 1 || budget <= 1 || !_rhs_flat_supported(env, dynamic_effectors)
        inner_ed = (active_sats <= 1 && budget > 1) ?
            effector_decision : _with_serial_effector_decision(effector_decision)
        return (
            mode=:satellite_batch,
            allotment=1,
            scheduler=:static,
            dominant_axis=:satellite,
            policy_applied=inner_ed.policy_applied,
            effector_decision=inner_ed,
        )
    end

    # Flat queue requires a minimum thread budget to amortise channel/worker overhead.
    if budget < env.flat_min_thread_budget
        return (
            mode=:satellite_batch,
            allotment=1,
            scheduler=:static,
            dominant_axis=:satellite,
            policy_applied=true,
            effector_decision=_with_serial_effector_decision(effector_decision),
        )
    end

    shared_buffers = _effector_shared_buffers(p)
    per_effector_cost_ns = _effector_observed_cost_ns_per_item(env, shared_buffers)
    estimated_total_work_ns   = per_effector_cost_ns * active_sats * n_effectors
    estimated_work_per_worker = estimated_total_work_ns / max(1, budget)
    # Flat queue is profitable only when total work is large AND each worker gets
    # enough to pay for dispatch/wakeup overhead (§5 per-worker gate).
    many_heavy_effectors = _rhs_effectors_have_heavy_or_heterogeneous_cost(dynamic_effectors, env.flat_cost_heterogeneity_threshold) &&
        estimated_total_work_ns   >= env.flat_work_ns_threshold &&
        estimated_work_per_worker >= env.flat_work_per_worker_ns_threshold

    if active_sats >= env.flat_min_sats &&
       (n_effectors >= env.flat_min_effectors || _rhs_flat_has_batch_privileged_effector(dynamic_effectors)) &&
       many_heavy_effectors
        # Same nested-outer-split hazard as the harmonics-batch/single-invsq
        # routes above.
        if outer_active && !env.harmonics_batch_allow_with_outer
            return (
                mode=:satellite_batch,
                allotment=1,
                scheduler=:static,
                dominant_axis=:satellite,
                policy_applied=true,
                effector_decision=_with_serial_effector_decision(effector_decision),
            )
        end
        return (
            mode=:flat_constellation_effector_queue,
            allotment=min(max(1, budget), active_sats * n_effectors),
            scheduler=:dynamic,
            dominant_axis=:flat_effector,
            policy_applied=true,
            effector_decision=_with_serial_effector_decision(effector_decision),
        )
    end

    # Satellite batch saturates the thread pool: every available thread handles at
    # least one satellite, so nested per-effector threads add only overhead.
    if _satellite_batch_saturates_pool(active_sats, budget)
        return (
            mode=:satellite_batch,
            allotment=1,
            scheduler=:static,
            dominant_axis=:satellite,
            policy_applied=true,
            effector_decision=_with_serial_effector_decision(effector_decision),
        )
    end

    # Few satellites, many threads: each satellite gets its own effector-reduce.
    return (
        mode=:per_satellite_effector_reduce,
        allotment=effector_decision.allotment,
        scheduler=:auto,
        dominant_axis=:per_satellite_inner_effector,
        policy_applied=effector_decision.policy_applied,
        effector_decision=effector_decision,
    )
end

# drag_cache/lift_cache/cross_cache in SaveCache all start as empty Vectors and grow
# lazily via _store_vector_cache!'s `resize!` (aerodynamic_wrench_models.jl). That's
# only safe called sequentially: AerodynamicCoefficientfM's wrench is dispatched
# across persistent worker threads by the flat-constellation-effector-queue route,
# so two satellites landing on different workers can both see the cache too short
# and both call `resize!` on the same shared Vector concurrently. Julia 1.12 catches
# this and throws ConcurrencyViolationError at small satellite counts; at larger
# counts (tighter timing races) it has corrupted memory badly enough to segfault
# instead. Pre-sizing to num_sats here, before any threaded dispatch can begin,
# means every per-satellite write lands on an already-appropriately-sized Vector at
# a distinct index, so the racy resize path is never reached.
function _initialize_save_cache_buffers!(p)::Nothing
    n_sats = length(p.args.dynamics_model.spacecraft)
    zero_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    for cache in (p.save_cache.drag_cache, p.save_cache.lift_cache, p.save_cache.cross_cache)
        if length(cache) != n_sats
            resize!(cache, n_sats)
        end
        fill!(cache, zero_vec)
    end
    return nothing
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

# in_atmosphere[] otherwise defaults to false for every satellite (runtime_types.jl)
# and is only ever flipped by the up/down-crossing event callback in
# event_callbacks.jl. A satellite whose initial orbit never crosses EI --
# because it starts (and stays) below it, e.g. a circular low-altitude orbit --
# would then incorrectly read as "not in atmosphere" for the entire mission,
# silently skipping the vacuum-predicted GRAM cache and the finer
# dt_max_atmosphere step size. Set the flag from the actual starting altitude
# instead of leaving every satellite to default to the "above the atmosphere"
# state regardless of where it actually starts.
function _initialize_in_atmosphere_flags!(p, initial_conditions)::Nothing
    sc_state = initial_conditions.sc
    n = length(sc_state)
    length(p.shared_buffers.in_atmosphere) == n || resize!(p.shared_buffers.in_atmosphere, n)
    planet = p.args.environment_model.planet
    ei_m = p.args.environment_model.EI * 1e3
    @inbounds for i in 1:n
        alt = norm(_state_position_ii(initial_conditions, i)) - planet.Rp_e
        p.shared_buffers.in_atmosphere[i] = alt <= ei_m
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
        empty!(memo.body_positions_j2000_m)
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
        if _is_nbody_effector_like(effector) && !isempty(effector.body_names)
            return true
        end
    end
    return false
end

function _collect_nbody_query_names(dynamic_effectors::Tuple)::Vector{String}
    names = String[]
    seen = Set{String}()
    @inbounds for effector in dynamic_effectors
        if !_is_nbody_effector_like(effector)
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

@inline function _ephemerides_time_seconds_flexible(initial_time, ephemerides_model)::Float64
    if applicable(SimulationModel.ephemerides_time_seconds, initial_time, ephemerides_model)
        return SimulationModel.ephemerides_time_seconds(initial_time, ephemerides_model)
    end

    model_name = nameof(typeof(ephemerides_model))
    if model_name == :SpiceEphemeridesModel
        start_epoch = SimulationModel.EphemeridesModels.from_utc(
            SimulationModel.EphemeridesModels._initial_time_datetime(initial_time)
        )
        return lock(RuntimeServices.SPICE_LOCK) do
            utc2et(SimulationModel.EphemeridesModels.to_utc(start_epoch))
        end
    elseif model_name == :SimpleEphemeridesModel
        start_time = SimulationModel.EphemeridesModels._initial_time_datetime(initial_time)
        elapsed_ms = Dates.value(start_time - SimulationModel.EphemeridesModels._J2000_UTC)
        return Float64(ephemerides_model.reference_epoch_seconds) + elapsed_ms / 1000.0
    end

    throw(MethodError(SimulationModel.ephemerides_time_seconds, (initial_time, ephemerides_model)))
end

@inline function _nbody_ephemeris_cache_sample_count(mission_end_s::Float64, dt_s::Float64)::Int
    return max(2, Int(ceil(mission_end_s / dt_s)) + 1)
end

function _nbody_ephemeris_body_index_by_name(body_query_names::Vector{String})::Dict{String, Int}
    body_index_by_name = Dict{String, Int}()
    @inbounds for (idx, body_name) in pairs(body_query_names)
        body_index_by_name[body_name] = idx
    end
    return body_index_by_name
end

function _nbody_ephemeris_cache_from_samples(
    primary_body_name::String,
    body_query_names::Vector{String},
    ets::Vector{Float64},
    positions_j2000_m::Matrix{SVector{3, Float64}}
)::SimulationModel.NBodyEphemerisCache
    body_index_by_name = _nbody_ephemeris_body_index_by_name(body_query_names)
    return SimulationModel.NBodyEphemerisCache(
        primary_body_name,
        body_query_names,
        body_index_by_name,
        ets,
        positions_j2000_m
    )
end

@inline function _increment_atomic_counter!(counter)
    counter === nothing && return nothing
    Base.Threads.atomic_add!(counter, 1)
    return nothing
end

function _build_nbody_ephemeris_cache(
    primary_body_name::String,
    body_query_names::Vector{String},
    et_start::Float64,
    mission_end_s::Float64,
    dt_s::Float64;
    counter=nothing
)::SimulationModel.NBodyEphemerisCache
    n_samples = _nbody_ephemeris_cache_sample_count(mission_end_s, dt_s)
    n_bodies = length(body_query_names)
    ets = Vector{Float64}(undef, n_samples)
    positions = Matrix{SVector{3, Float64}}(undef, n_samples, n_bodies)

    lock(RuntimeServices.SPICE_LOCK) do
        @inbounds for sample_idx in 1:n_samples
            et = et_start + min((sample_idx - 1) * dt_s, mission_end_s)
            ets[sample_idx] = et
            for body_idx in 1:n_bodies
                body_query = body_query_names[body_idx]
                positions[sample_idx, body_idx] = SimulationModel.EphemeridesModels._spice_position_j2000_m_unlocked(body_query, et, primary_body_name)
                _increment_atomic_counter!(counter)
            end
        end
    end

    return _nbody_ephemeris_cache_from_samples(primary_body_name, copy(body_query_names), ets, positions)
end

@inline function _prewarmed_nbody_ephemeris_lookup(primary_body_name::String, body_query_names::Vector{String}, et_start::Float64, mission_end_s::Float64, dt_s::Float64)
    reuse_key = _nbody_ephemeris_reuse_key(primary_body_name, body_query_names, et_start, mission_end_s, dt_s)
    return _ephemeris_reuse_lookup(_NBODY_EPHEMERIS_PREWARMED_CACHE, reuse_key)
end

@inline function _register_prewarmed_nbody_ephemeris_cache!(
    primary_body_name::String,
    body_query_names::Vector{String},
    et_start::Float64,
    mission_end_s::Float64,
    dt_s::Float64,
    cache::SimulationModel.NBodyEphemerisCache;
    replace::Bool=false
)::SimulationModel.NBodyEphemerisCache
    reuse_key = _nbody_ephemeris_reuse_key(primary_body_name, body_query_names, et_start, mission_end_s, dt_s)
    return _ephemeris_explicit_cache_store!(_NBODY_EPHEMERIS_PREWARMED_CACHE, reuse_key, cache; replace=replace)
end

function _nbody_ephemeris_cache_payload(
    cache::SimulationModel.NBodyEphemerisCache,
    et_start::Float64,
    mission_end_s::Float64,
    dt_s::Float64
)
    return (
        schema_version=NBODY_EPHEMERIS_CACHE_SCHEMA_VERSION,
        created_utc=string(now(UTC)),
        primary_body_name=cache.primary_body_name,
        body_query_names=copy(cache.body_query_names),
        et_start=Float64(et_start),
        mission_end_s=Float64(mission_end_s),
        dt_s=Float64(dt_s),
        ets=copy(cache.ets),
        positions_j2000_m=copy(cache.positions_j2000_m),
    )
end

function _write_nbody_ephemeris_cache_file!(
    path::String,
    cache::SimulationModel.NBodyEphemerisCache,
    et_start::Float64,
    mission_end_s::Float64,
    dt_s::Float64
)::String
    payload = _nbody_ephemeris_cache_payload(cache, et_start, mission_end_s, dt_s)
    _atomic_write_file(path, tmp -> open(tmp, "w") do io
        serialize(io, payload)
    end)
    return path
end

@inline function _payload_field(payload, field::Symbol)
    hasproperty(payload, field) || throw(ArgumentError("N-body ephemeris cache payload missing required field $(repr(field))."))
    return getproperty(payload, field)
end

function _cache_from_nbody_ephemeris_payload(payload)
    payload isa NamedTuple || throw(ArgumentError("N-body ephemeris cache payload must be a NamedTuple serialized by SpaceAGORA."))

    schema_version = String(_payload_field(payload, :schema_version))
    schema_version == NBODY_EPHEMERIS_CACHE_SCHEMA_VERSION || throw(ArgumentError(
        "Unsupported N-body ephemeris cache schema_version=$(repr(schema_version)); expected $(repr(NBODY_EPHEMERIS_CACHE_SCHEMA_VERSION))."
    ))

    primary_body_name = String(_payload_field(payload, :primary_body_name))
    body_query_names = String[String(name) for name in _payload_field(payload, :body_query_names)]
    et_start = Float64(_payload_field(payload, :et_start))
    mission_end_s = Float64(_payload_field(payload, :mission_end_s))
    dt_s = Float64(_payload_field(payload, :dt_s))
    ets = Vector{Float64}(_payload_field(payload, :ets))
    positions_j2000_m = Matrix{SVector{3, Float64}}(_payload_field(payload, :positions_j2000_m))

    isfinite(et_start) || throw(ArgumentError("Serialized N-body ephemeris cache et_start must be finite."))
    (isfinite(mission_end_s) && mission_end_s > 0.0) || throw(ArgumentError("Serialized N-body ephemeris cache mission_end_s must be > 0."))
    (isfinite(dt_s) && dt_s > 0.0) || throw(ArgumentError("Serialized N-body ephemeris cache dt_s must be > 0."))
    isempty(body_query_names) && throw(ArgumentError("Serialized N-body ephemeris cache must include at least one body_query_name."))
    length(ets) >= 2 || throw(ArgumentError("Serialized N-body ephemeris cache must include at least two samples."))
    size(positions_j2000_m, 1) == length(ets) || throw(ArgumentError("Serialized N-body ephemeris cache positions sample count does not match ets length."))
    size(positions_j2000_m, 2) == length(body_query_names) || throw(ArgumentError("Serialized N-body ephemeris cache body count does not match body_query_names length."))

    cache = _nbody_ephemeris_cache_from_samples(primary_body_name, body_query_names, ets, positions_j2000_m)
    return (
        cache=cache,
        et_start=et_start,
        mission_end_s=mission_end_s,
        dt_s=dt_s,
    )
end

function _load_nbody_ephemeris_cache!(path::String; replace::Bool=true)::SimulationModel.NBodyEphemerisCache
    payload = open(path, "r") do io
        deserialize(io)
    end
    loaded = _cache_from_nbody_ephemeris_payload(payload)
    cache = _register_prewarmed_nbody_ephemeris_cache!(
        loaded.cache.primary_body_name,
        loaded.cache.body_query_names,
        loaded.et_start,
        loaded.mission_end_s,
        loaded.dt_s,
        loaded.cache;
        replace=replace
    )
    return cache
end

function _prewarm_nbody_ephemeris_cache(
    args;
    dt_s::Union{Nothing, Real}=nothing,
    mission_end_s::Union{Nothing, Real}=nothing,
    save_path::Union{Nothing, AbstractString}=nothing
)::SimulationModel.NBodyEphemerisCache
    _validate_ephemerides_support!(args)
    _has_active_nbody_effector(args.dynamics_model.dynamic_effectors) || throw(ArgumentError(
        "prewarm_nbody_ephemeris_cache requires at least one active NBodyGravityModel with non-empty body_names."
    ))

    resolved_dt_s = isnothing(dt_s) ? _nbody_ephemeris_cache_dt_s() : Float64(dt_s)
    (isfinite(resolved_dt_s) && resolved_dt_s > 0.0) || throw(ArgumentError("dt_s must be > 0.0, got $(repr(dt_s))."))
    resolved_mission_end_s = isnothing(mission_end_s) ? args.mission_configuration.mission_time : Float64(mission_end_s)
    (isfinite(resolved_mission_end_s) && resolved_mission_end_s > 0.0) || throw(ArgumentError("mission_end_s must be > 0.0, got $(repr(mission_end_s))."))

    n_samples = _nbody_ephemeris_cache_sample_count(resolved_mission_end_s, resolved_dt_s)
    max_samples = _nbody_ephemeris_cache_max_samples()
    if n_samples > max_samples
        throw(ArgumentError(
            "N-body ephemeris prewarm requires samples=$n_samples, which exceeds SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES=$max_samples. " *
            "Increase the max-samples policy or choose a larger dt_s / shorter mission_end_s."
        ))
    end

    body_query_names = _collect_nbody_query_names(args.dynamics_model.dynamic_effectors)
    isempty(body_query_names) && throw(ArgumentError(
        "prewarm_nbody_ephemeris_cache requires at least one third body to cache."
    ))

    ephemerides_model = args.environment_model.ephemerides_model
    et_start = _ephemerides_time_seconds_flexible(args.initial_time, ephemerides_model)
    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(args.environment_model.planet.name)

    cache = _prewarmed_nbody_ephemeris_lookup(primary_body_name, body_query_names, et_start, resolved_mission_end_s, resolved_dt_s)
    if !(cache isa SimulationModel.NBodyEphemerisCache)
        reused = _ephemeris_reuse_lookup(
            _NBODY_EPHEMERIS_REUSE_CACHE,
            _nbody_ephemeris_reuse_key(primary_body_name, body_query_names, et_start, resolved_mission_end_s, resolved_dt_s)
        )
        if reused isa SimulationModel.NBodyEphemerisCache
            cache = _register_prewarmed_nbody_ephemeris_cache!(
                primary_body_name,
                body_query_names,
                et_start,
                resolved_mission_end_s,
                resolved_dt_s,
                reused
            )
        else
            cache = _build_nbody_ephemeris_cache(primary_body_name, body_query_names, et_start, resolved_mission_end_s, resolved_dt_s)
            cache = _register_prewarmed_nbody_ephemeris_cache!(
                primary_body_name,
                body_query_names,
                et_start,
                resolved_mission_end_s,
                resolved_dt_s,
                cache
            )
        end
    end

    if save_path !== nothing
        path_str = String(save_path)
        isempty(strip(path_str)) && throw(ArgumentError("save_path must not be empty when provided."))
        _write_nbody_ephemeris_cache_file!(path_str, cache, et_start, resolved_mission_end_s, resolved_dt_s)
    end
    return cache
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

    lock(RuntimeServices.SPICE_LOCK) do
        @inbounds for sample_idx in 1:n_samples
            et = et_start + min((sample_idx - 1) * dt_s, mission_end_s)
            ets[sample_idx] = et
            positions[sample_idx] = SimulationModel.EphemeridesModels._spice_position_j2000_m_unlocked("sun", et, primary_body_name)
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
    n_samples = _nbody_ephemeris_cache_sample_count(mission_end_s, dt_s)
    max_samples = _nbody_ephemeris_cache_max_samples()
    if n_samples > max_samples
        @warn "N-body ephemeris cache disabled: required samples=$n_samples exceeds SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES=$max_samples."
        return nothing
    end

    primary_body_name = SimulationModel.DynamicEffectors._spice_query_name(p.args.environment_model.planet.name)
    prewarmed = _prewarmed_nbody_ephemeris_lookup(primary_body_name, body_query_names, et_start, mission_end_s, dt_s)
    if prewarmed isa SimulationModel.NBodyEphemerisCache
        p.shared_buffers.nbody_ephemeris_cache[] = prewarmed
        return nothing
    end
    if _ephemeris_reuse_enabled()
        reuse_key = _nbody_ephemeris_reuse_key(primary_body_name, body_query_names, et_start, mission_end_s, dt_s)
        reused = _ephemeris_reuse_lookup(_NBODY_EPHEMERIS_REUSE_CACHE, reuse_key)
        if reused isa SimulationModel.NBodyEphemerisCache
            p.shared_buffers.nbody_ephemeris_cache[] = reused
            return nothing
        end
    end
    n_bodies = length(body_query_names)
    cache_value = _build_nbody_ephemeris_cache(
        primary_body_name,
        body_query_names,
        et_start,
        mission_end_s,
        dt_s;
        counter=p.shared_buffers.spice_runtime_counters.nbody_spkpos_cache_build_calls
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
