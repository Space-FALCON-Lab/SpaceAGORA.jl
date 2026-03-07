@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline _gram_track_cache_ignore_time_window() = _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW", true)
@inline _gram_track_cache_target_use_j2() = _parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_TARGET_USE_J2", true)

@inline function _gram_entry_target_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_ENTRY_TARGET_MODE", "allen_eggers")))
    if mode in ("off", "none", "0", "false", "no")
        return :off
    elseif mode in ("allen_eggers", "allen-eggers", "allen", "ae", "on", "1", "true", "yes", "auto")
        return :allen_eggers
    end
    throw(ArgumentError("Unsupported SPACEAGORA_GRAM_ENTRY_TARGET_MODE='$mode'. Use one of: off, allen_eggers."))
end

@inline _gram_entry_target_cd() = max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_CD", 1.5))
@inline _gram_entry_target_dt() = max(0.05, _parse_float_env("SPACEAGORA_GRAM_ENTRY_TARGET_DT_S", 0.5))
@inline function _gram_entry_target_max_steps()::Int
    raw = strip(get(ENV, "SPACEAGORA_GRAM_ENTRY_TARGET_MAX_STEPS", "400"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_GRAM_ENTRY_TARGET_MAX_STEPS must be an integer value, got '$raw'"))
    end
    return max(8, parsed)
end

@inline function _density_callback_parallel_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL")
end

@inline function _density_callback_thread_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD", 8)
end

@inline function _density_callback_allow_with_outer()::Bool
    return _parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)
end

@inline function _density_batch_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_DENSITY_BATCH_PARALLEL")
end

@inline function _density_batch_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_DENSITY_BATCH_THRESHOLD", 2)
end

@inline function _density_batch_enabled(num_sats::Int)::Bool
    mode = _density_batch_mode()
    if mode == :off
        return false
    elseif mode == :on
        return num_sats > 0
    end
    return num_sats >= _density_batch_threshold()
end

@inline function _gram_isolated_pool_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_GRAM_ISOLATED_POOL"; default="off")
end

@inline function _gram_isolated_pool_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD", 4)
end

@inline function _gram_isolated_pool_max_workers()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS", max(1, Threads.nthreads()))
end

@inline function _gram_isolated_pool_enabled(num_items::Int)::Bool
    mode = _gram_isolated_pool_mode()
    if mode == :off
        return false
    elseif mode == :on
        return num_items > 0
    end
    return Threads.nthreads() > 1 && num_items >= _gram_isolated_pool_threshold()
end

@inline function _control_callback_parallel_mode()::Symbol
    return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL")
end

@inline function _control_callback_thread_threshold()::Int
    return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD", 8)
end

@inline function _control_callback_allow_with_outer()::Bool
    return _parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)
end

@inline function _thermal_callback_parallel_mode()::Symbol
    if haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_PARALLEL")
        return ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_THERMAL_CALLBACK_PARALLEL")
    end
    return _density_callback_parallel_mode()
end

@inline function _thermal_callback_thread_threshold()::Int
    if haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD")
        return ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD", 8)
    end
    return _density_callback_thread_threshold()
end

@inline function _thermal_callback_allow_with_outer()::Bool
    if haskey(ENV, "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER")
        return _parse_bool_env("SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)
    end
    return _density_callback_allow_with_outer()
end

@inline function _callback_outer_parallel_hint()::Bool
    return ParallelPolicy.outer_parallel_active()
end

# Extend this for custom user density models as needed:
# SimulationModel.SimulationCallbacks.density_model_threadsafe(::MyDensityModel) = true
@inline density_model_threadsafe(::AbstractDensityModel)::Bool = false
@inline density_model_threadsafe(::NoAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.ExponentialAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.PolynomialFitAtmosphereModel)::Bool = true
# GRAM C-wrapper calls are serialized inside getDensity via SimulationModel.GRAM_LOCK.
@inline density_model_threadsafe(::EnvironmentModels.GRAMAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.GRAMAtmosphereModelSurrogate)::Bool = true

@inline _is_gram_density_model(model)::Bool =
    model isa EnvironmentModels.GRAMAtmosphereModel ||
    model isa EnvironmentModels.GRAMAtmosphereModelSurrogate

@inline function _gram_track_trajectory_supported(density_model)::Bool
    _is_gram_density_model(density_model) || return false
    hasproperty(density_model, :gram) || return false
    hasproperty(density_model, :gram_atmosphere) || return false
    gram_driver = try
        getproperty(density_model, :gram)
    catch
        return false
    end
    return hasproperty(gram_driver, :generate_trajectory)
end

@inline function _density_callback_thread_decision(args::SimulationConfiguration, num_sats::Int)
    mode = _density_callback_parallel_mode()
    outer_active = _callback_outer_parallel_hint()
    allow_with_outer = _density_callback_allow_with_outer()

    model = args.environment_model.density_model
    model_threadsafe = density_model_threadsafe(model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_ASSUME_THREADSAFE", false)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=_density_callback_thread_threshold(),
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:density_callback
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end

@inline function _density_callback_use_threads(args::SimulationConfiguration, num_sats::Int)::Bool
    return _density_callback_thread_decision(args, num_sats).use_threads
end

# Extend this for custom user control models as needed:
# SimulationModel.SimulationCallbacks.control_model_threadsafe(::MyControlModel) = true
@inline control_model_threadsafe(::Any)::Bool = false
@inline control_model_threadsafe(::BaseThrusterModel)::Bool = true

@inline function _control_callback_thread_decision(control_model, num_sats::Int, use_invokelatest::Bool)
    if use_invokelatest
        return (use_threads=false, allotment=1, mode=:off, policy_applied=false)
    end
    mode = _control_callback_parallel_mode()
    outer_active = _callback_outer_parallel_hint()
    allow_with_outer = _control_callback_allow_with_outer()

    model_threadsafe = control_model_threadsafe(control_model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE", false)
        return (use_threads=false, allotment=1, mode=mode, policy_applied=false)
    end
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=_control_callback_thread_threshold(),
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:control_callback
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end

@inline function _control_callback_use_threads(control_model, num_sats::Int, use_invokelatest::Bool)::Bool
    return _control_callback_thread_decision(control_model, num_sats, use_invokelatest).use_threads
end

@inline function _thermal_callback_thread_decision(num_sats::Int)
    mode = _thermal_callback_parallel_mode()
    outer_active = _callback_outer_parallel_hint()
    allow_with_outer = _thermal_callback_allow_with_outer()
    policy = ParallelPolicy.thread_policy_decision(
        num_sats;
        mode=mode,
        threshold=_thermal_callback_thread_threshold(),
        outer_active=outer_active,
        allow_with_outer=allow_with_outer,
        source=:thermal_callback
    )
    return (use_threads=policy.use_threads, allotment=policy.allotment, mode=mode, policy_applied=true)
end

@inline function _density_model_for_sat(
    density_models::AbstractVector{<:AbstractDensityModel},
    fallback_model::AbstractDensityModel,
    sat_idx::Int
)
    if sat_idx <= length(density_models)
        return density_models[sat_idx]
    end
    return fallback_model
end

@inline function _density_model_for_sat(p, sat_idx::Int)
    return _density_model_for_sat(
        p.shared_buffers.density_models,
        p.args.environment_model.density_model,
        sat_idx,
    )
end

@inline function _density_batch_model_for_callback(
    density_models::AbstractVector{<:AbstractDensityModel},
    fallback_model::AbstractDensityModel,
    num_sats::Int
)
    if isempty(density_models)
        return fallback_model
    end
    if length(density_models) < num_sats
        return nothing
    end
    first_model = density_models[1]
    @inbounds for i in 2:num_sats
        density_models[i] === first_model || return nothing
    end
    return first_model
end

@inline function _density_batch_model_for_callback(p, num_sats::Int)
    return _density_batch_model_for_callback(
        p.shared_buffers.density_models,
        p.args.environment_model.density_model,
        num_sats,
    )
end

@inline function _gram_isolated_pool_batch_model_for_callback(
    density_models::AbstractVector{<:AbstractDensityModel},
    fallback_model::AbstractDensityModel,
    num_sats::Int
)
    _gram_isolated_pool_enabled(num_sats) || return nothing
    isempty(density_models) || return nothing
    model = fallback_model
    return model isa EnvironmentModels.GRAMAtmosphereModel ? model : nothing
end

@inline function _gram_isolated_pool_batch_model_for_callback(p, num_sats::Int)
    return _gram_isolated_pool_batch_model_for_callback(
        p.shared_buffers.density_models,
        p.args.environment_model.density_model,
        num_sats,
    )
end

@inline function _gram_batch_elapsed_time(el_time::Float64, idx::Int)::Float64
    return el_time
end

@inline function _gram_batch_elapsed_time(el_time::AbstractVector{<:Real}, idx::Int)::Float64
    return Float64(el_time[idx])
end

@inline function _gram_isolated_pool_density_state(
    model::EnvironmentModels.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p,
    model_lock::ReentrantLock
)::Tuple{Float64, Float64, SVector{3, Float64}}
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0
    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return EnvironmentModels.density_polyfit(h, p)
    end
    return GRAMSuite.density_state(
        model.core,
        h,
        lat,
        lon,
        el_time,
        wind;
        lock_obj=model_lock,
        vacuum_temperature=p.args.environment_model.planet.T_ref
    )
end

@inline function _gram_density_cache_for_sat!(
    caches::Vector{Union{Nothing, GramTrackCache}},
    sat_idx::Int
)::GramTrackCache
    if sat_idx <= length(caches)
        cache = @inbounds caches[sat_idx]
        if cache === nothing
            cache = GramTrackCache()
            @inbounds caches[sat_idx] = cache
        end
        return cache
    end
    return GramTrackCache()
end

function _ensure_gram_isolated_pool!(
    p,
    template_model::EnvironmentModels.GRAMAtmosphereModel,
    workers::Int
)::Tuple{Vector{EnvironmentModels.GRAMAtmosphereModel}, Vector{ReentrantLock}}
    shared = p.shared_buffers
    models = shared.gram_isolated_pool_models
    locks = shared.gram_isolated_pool_locks
    if workers <= 0
        return models, locks
    end
    rebuild = length(models) != workers || length(locks) != workers
    if rebuild
        empty!(models)
        empty!(locks)
        sizehint!(models, workers)
        sizehint!(locks, workers)
        @inbounds for _ in 1:workers
            push!(models, deepcopy(template_model))
            push!(locks, ReentrantLock())
        end
    end
    return models, locks
end

@inline function _gram_isolated_pool_batch_eval!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    density_model,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p;
    allotment_hint::Int=max(1, Threads.nthreads())
)::Bool
    return false
end

function _gram_isolated_pool_batch_eval!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    density_model::EnvironmentModels.GRAMAtmosphereModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p;
    allotment_hint::Int=max(1, Threads.nthreads())
)::Bool
    n = length(hs)
    _gram_isolated_pool_enabled(n) || return false
    length(rhos) == n || return false
    length(Ts) == n || return false
    length(winds) == n || return false
    length(lats) == n || return false
    length(lons) == n || return false
    if el_time isa AbstractVector{<:Real}
        length(el_time) == n || return false
    end

    max_allotment = min(max(1, allotment_hint), _gram_isolated_pool_max_workers())
    workers = ParallelPolicy.thread_worker_count(n, max_allotment)
    workers > 1 || return false

    models, locks = _ensure_gram_isolated_pool!(p, density_model, workers)
    ParallelPolicy.threaded_foreach_worker(n, workers) do worker_id, idx
        model_i = models[worker_id]::EnvironmentModels.GRAMAtmosphereModel
        lock_i = locks[worker_id]
        h = Float64(hs[idx])
        lat = Float64(lats[idx])
        lon = Float64(lons[idx])
        t_i = _gram_batch_elapsed_time(el_time, idx)
        rho_i, T_i, wind_i = _gram_isolated_pool_density_state(model_i, h, lat, lon, t_i, wind, p, lock_i)
        @inbounds begin
            rhos[idx] = rho_i
            Ts[idx] = T_i
            winds[idx] = wind_i
        end
    end
    return true
end


@inline function _uses_atmospheric_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
            return true
        end
    end
    return false
end

@inline function _uses_j2_gravity_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa InverseSquaredJ2GravityModel
            return true
        end
    end
    return false
end

@inline function _requires_density_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _uses_atmospheric_dynamic_effector(effectors) || !(args.environment_model.density_model isa NoAtmosphereModel)
end

@inline function _requires_orbit_end_callback(args::SimulationConfiguration)::Bool
    return args.mission_configuration.mission_type == MissionOrbits
end

@inline function _entry_target_count()::Int
    raw = strip(get(ENV, "SPACEAGORA_ENTRY_TARGET_COUNT", "0"))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_ENTRY_TARGET_COUNT must be an integer value, got '$raw'"))
    end
    parsed >= 0 || throw(ArgumentError("SPACEAGORA_ENTRY_TARGET_COUNT must be >= 0, got $parsed"))
    return parsed
end

@inline function _requires_entry_end_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    return _entry_target_count() > 0 && _requires_density_callback(effectors, args)
end

@inline function _requires_drag_state_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    if !_requires_density_callback(effectors, args)
        return false
    end
    tol = args.integration_tolerances
    return tol.dt_max_atmosphere != tol.dt_max_orbit ||
           tol.reltol_atmosphere != tol.reltol_orbit ||
           tol.abstol_atmosphere != tol.abstol_orbit
end

@inline function _requires_thermal_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    # Thermal-rate evaluation requires atmospheric state (ρ, T, wind), populated by density callback.
    return _requires_density_callback(effectors, args)
end

@inline function _requires_quaternion_projection_callback(args::SimulationConfiguration)::Bool
    return args.mission_configuration.orientation_sim
end

@inline function _resolved_component_tolerance(component_tol::Float64, baseline_tol::Float64)::Float64
    return component_tol == 0.0 ? baseline_tol : component_tol
end

function _callback_tolerances_for_phase(template_reltol, template_abstol, args::SimulationConfiguration, in_atmosphere::Bool)
    tol = args.integration_tolerances
    baseline_reltol = in_atmosphere ? tol.reltol_atmosphere : tol.reltol_orbit
    baseline_abstol = in_atmosphere ? tol.abstol_atmosphere : tol.abstol_orbit
    if template_reltol isa Number && template_abstol isa Number
        return baseline_reltol, baseline_abstol
    end

    reltol_mass = _resolved_component_tolerance(tol.reltol_mass, baseline_reltol)
    abstol_mass = _resolved_component_tolerance(tol.abstol_mass, baseline_abstol)
    reltol_heat = _resolved_component_tolerance(tol.reltol_heat_load, baseline_reltol)
    abstol_heat = _resolved_component_tolerance(tol.abstol_heat_load, baseline_abstol)
    reltol_ω = _resolved_component_tolerance(tol.reltol_angular_rate, baseline_reltol)
    abstol_ω = _resolved_component_tolerance(tol.abstol_angular_rate, baseline_abstol)

    reltol_new = copy(template_reltol)
    abstol_new = copy(template_abstol)
    reltol_new .= baseline_reltol
    abstol_new .= baseline_abstol

    @inbounds for i in eachindex(reltol_new.sc)
        reltol_new.sc[i].mass = reltol_mass
        abstol_new.sc[i].mass = abstol_mass
        reltol_new.sc[i].heat_loads .= reltol_heat
        abstol_new.sc[i].heat_loads .= abstol_heat
        if hasproperty(reltol_new.sc[i], :q)
            reltol_new.sc[i].q .= tol.reltol_quaternion
            abstol_new.sc[i].q .= tol.abstol_quaternion
        end
        if hasproperty(reltol_new.sc[i], :ω)
            reltol_new.sc[i].ω .= reltol_ω
            abstol_new.sc[i].ω .= abstol_ω
        end
    end
    return reltol_new, abstol_new
end

function get_callbacks(
    num_sats::Int,
    effectors::Tuple,
    args::SimulationConfiguration;
    saved_values=nothing,
    save_fields=nothing
)::CallbackSet
    save_fields_resolved = _resolve_save_fields(save_fields, args)
    callbacks = Any[
        get_impact_callback(num_sats),
        update_planet_frame_callback()
    ]

    if _requires_density_callback(effectors, args)
        push!(callbacks, get_density_callback(num_sats, effectors, args))
    end

    if _requires_thermal_callback(effectors, args)
        push!(callbacks, get_thermal_callback(num_sats, args))
    end

    if _requires_orbit_end_callback(args)
        push!(callbacks, get_orbit_end_callback(num_sats))
    end

    if _requires_entry_end_callback(effectors, args)
        push!(callbacks, get_entry_end_callback(num_sats, args))
    end

    if _requires_drag_state_callback(effectors, args)
        push!(callbacks, get_drag_state_callback(num_sats))
    end

    append!(callbacks, get_navigation_callbacks(num_sats, args))
    append!(callbacks, get_control_callbacks(num_sats, args))
    append!(callbacks, get_guidance_callbacks(num_sats, args))
    if _requires_quaternion_projection_callback(args)
        push!(callbacks, get_quaternion_projection_callback(num_sats, args))
    end
    push!(callbacks, get_data_saving_callback(num_sats, args, save_fields_resolved, saved_values))

    return CallbackSet(callbacks...)
end

@inline function _planet_lpi_from_backend(planet, ephemerides_model, et::Float64, counter::Base.Threads.Atomic{Int64})::SMatrix{3, 3, Float64}
    if ephemerides_requires_spice(ephemerides_model)
        Base.Threads.atomic_add!(counter, 1)
    end
    return planet_frame_lpi(planet, et, ephemerides_model)
end

@inline function _planet_lpi_from_cache(cache::PlanetFrameEphemerisCache, et::Float64)::Union{Nothing, SMatrix{3, 3, Float64}}
    ets = cache.ets
    n_samples = length(ets)
    n_samples >= 2 || return nothing
    if et < ets[1] || et > ets[n_samples]
        return nothing
    end

    idx = searchsortedlast(ets, et)
    if idx <= 0
        return nothing
    elseif idx >= n_samples
        return rot(cache.quaternions[n_samples])
    end

    et0 = ets[idx]
    et1 = ets[idx + 1]
    if et1 <= et0
        return rot(cache.quaternions[idx])
    end

    α = (et - et0) / (et1 - et0)
    q0 = cache.quaternions[idx]
    q1 = cache.quaternions[idx + 1]
    # Keep interpolation on one quaternion hemisphere to avoid discontinuities.
    if dot(q0, q1) < 0.0
        q1 = -q1
    end
    q_interp = normalize((1.0 - α) * q0 + α * q1)
    return rot(q_interp)
end

function update_planet_frame_callback()
    condition(u, t, integrator) = true
    et_start = Ref{Float64}(0.0) # Store the start epoch in a Ref so it can be updated in the affect! function
    function affect!(integrator)
        p = integrator.p
        planet = p.args.environment_model.planet
        ephemerides_model = p.args.environment_model.ephemerides_model
        et = et_start[] + integrator.t
        pxform_counter = p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls
        cache_entry = p.shared_buffers.planet_frame_ephemeris_cache[]
        l_pi = if cache_entry isa PlanetFrameEphemerisCache
            cached = _planet_lpi_from_cache(cache_entry, et)
            cached === nothing ? _planet_lpi_from_backend(planet, ephemerides_model, et, pxform_counter) : cached
        else
            _planet_lpi_from_backend(planet, ephemerides_model, et, pxform_counter)
        end
        planet.L_PI .= l_pi
    end

    function init_affect!(cb, u, t, integrator)
        p = integrator.p
        ephemerides_model = p.args.environment_model.ephemerides_model
        et_start[] = ephemerides_time_seconds(p.args.initial_time, ephemerides_model)
        affect!(integrator) # Call the affect! function at the start of the simulation to initialize the planet frame
    end

    return DiscreteCallback(condition, affect!, initialize=init_affect!)
end
# Factory function to build the callback
function get_density_callback(num_sats::Int, args::SimulationConfiguration)
    return get_density_callback(num_sats, args.dynamics_model.dynamic_effectors, args)
end

function get_density_callback(num_sats::Int, effectors::Tuple, args::SimulationConfiguration)
    cache_cfg = _gram_track_cache_config()
    stats_enabled = _gram_runtime_stats_enabled()
    target_include_j2 = _gram_track_cache_target_use_j2() && _uses_j2_gravity_effector(effectors)
    function update_density_sat!(i::Int, p, u, t, segment_end_t::Float64, density_model, caches::Vector{Union{Nothing, GramTrackCache}})
        if stats_enabled
            _gram_runtime_stats_update!(s -> begin
                s.density_calls += 1
            end)
        end
        pos = SVector{3, Float64}(u.sc[i].pos)
        vel = SVector{3, Float64}(u.sc[i].vel)
        rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet)
        alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet)
        ρ, T, wind_vec = if _gram_track_cache_enabled(cache_cfg, density_model)
            if stats_enabled
                _gram_runtime_stats_update!(s -> begin
                    s.cache_enabled_calls += 1
                end)
            end
            cache_horizon_s, cache_alt_tol_m, cache_ang_tol_rad, cache_points = _gram_track_cache_profile(cache_cfg, p, alt)
            cache = _gram_density_cache_for_sat!(caches, i)
            seg = if stats_enabled
                seg0 = _gram_track_cache_segment(cache, t)
                if seg0 === nothing
                    _gram_runtime_stats_update!(s -> begin
                        s.cache_misses += 1
                        s.miss_time_window += 1
                    end)
                    nothing
                else
                    idx0, x0 = seg0
                    alt_interp = _lerp(cache.alts[idx0], cache.alts[idx0 + 1], x0)
                    lat_interp = _lerp_angle_rad(cache.lats[idx0], cache.lats[idx0 + 1], x0)
                    lon_interp = _lerp_angle_rad(cache.lons[idx0], cache.lons[idx0 + 1], x0)
                    alt_err = abs(alt - alt_interp)
                    lat_err = _angdiff_rad(lat, lat_interp)
                    lon_err = _angdiff_rad(lon, lon_interp)
                    _gram_runtime_stats_update!(s -> begin
                        s.state_error_samples += 1
                        s.alt_err_abs_max_m = max(s.alt_err_abs_max_m, alt_err)
                        s.lat_err_abs_max_deg = max(s.lat_err_abs_max_deg, rad2deg(lat_err))
                        s.lon_err_abs_max_deg = max(s.lon_err_abs_max_deg, rad2deg(lon_err))
                        s.alt_err_abs_sum_m += alt_err
                        s.lat_err_abs_sum_deg += rad2deg(lat_err)
                        s.lon_err_abs_sum_deg += rad2deg(lon_err)
                    end)
                    if abs(alt - alt_interp) <= cache_alt_tol_m &&
                       _angdiff_rad(lat, lat_interp) <= cache_ang_tol_rad &&
                       _angdiff_rad(lon, lon_interp) <= cache_ang_tol_rad
                        _gram_runtime_stats_update!(s -> begin
                            s.cache_hits += 1
                        end)
                        seg0
                    else
                        _gram_runtime_stats_update!(s -> begin
                            s.cache_misses += 1
                            s.miss_state_tolerance += 1
                        end)
                        nothing
                    end
                end
            else
                _gram_track_cache_ready(cache, t, alt, lat, lon, cache_alt_tol_m, cache_ang_tol_rad)
            end
            if seg !== nothing
                idx, x = seg
                _gram_track_cache_eval(cache, idx, x)
            else
                _gram_track_cache_refresh!(
                    cache,
                    density_model,
                    p,
                    pos,
                    vel,
                    alt,
                    lat,
                    lon,
                    t,
                    cache_horizon_s,
                    cache_points,
                    cache_alt_tol_m,
                    cache_ang_tol_rad,
                    cache_cfg.transition_band_m,
                    segment_end_t,
                    target_include_j2,
                    i,
                    Float64(u.sc[i].mass)
                )
            end
        else
            if stats_enabled
                _gram_runtime_stats_update!(s -> begin
                    s.direct_calls += 1
                end)
            end
            getDensity(density_model, alt, lat, lon, t, true, p)
        end
        p.shared_buffers.densities[i] = ρ
        p.shared_buffers.temperatures[i] = T
        p.shared_buffers.winds[i] = wind_vec
        return nothing
    end
    
    # Logic for when the callback should trigger
    # Returning true means it runs at every integration step
    condition(u, t, integrator) = true

    # Logic for what happens when it triggers
    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        density_models = p.shared_buffers.density_models
        fallback_density_model = p.args.environment_model.density_model
        gram_density_caches = p.shared_buffers.gram_density_cache
        decision = _density_callback_thread_decision(args, num_sats)
        use_threads = decision.use_threads
        use_batch = false
        batch_model = nothing
        use_gram_isolated_pool = false
        if _density_batch_enabled(num_sats)
            batch_model = _density_batch_model_for_callback(density_models, fallback_density_model, num_sats)
            if !(batch_model === nothing) && !_gram_track_cache_enabled(cache_cfg, batch_model)
                use_batch = true
            end
        end
        if !use_batch
            gram_pool_model = _gram_isolated_pool_batch_model_for_callback(density_models, fallback_density_model, num_sats)
            if !(gram_pool_model === nothing) && !_gram_track_cache_enabled(cache_cfg, gram_pool_model)
                use_batch = true
                batch_model = gram_pool_model
                use_gram_isolated_pool = true
            end
        elseif batch_model isa EnvironmentModels.GRAMAtmosphereModel
            use_gram_isolated_pool = true
        end
        started_ns = time_ns()
        segment_end_t = try
            Float64(integrator.sol.prob.tspan[2])
        catch
            integrator.t + cache_cfg.orbit_horizon_s
        end

        if use_batch
            if stats_enabled
                _gram_runtime_stats_update!(s -> begin
                    s.density_calls += num_sats
                    s.direct_calls += num_sats
                end)
            end
            alts = p.shared_buffers.density_batch_altitudes
            lats = p.shared_buffers.density_batch_latitudes
            lons = p.shared_buffers.density_batch_longitudes
            if use_threads
                ParallelPolicy.threaded_foreach_persistent(:density_callback, num_sats, decision.allotment) do i
                    @inbounds begin
                        pos = SVector{3, Float64}(u.sc[i].pos)
                        vel = SVector{3, Float64}(u.sc[i].vel)
                        rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet)
                        alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet)
                        alts[i] = alt
                        lats[i] = lat
                        lons[i] = lon
                    end
                end
            else
                @inbounds for i in 1:num_sats
                    pos = SVector{3, Float64}(u.sc[i].pos)
                    vel = SVector{3, Float64}(u.sc[i].vel)
                    rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet)
                    alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet)
                    alts[i] = alt
                    lats[i] = lat
                    lons[i] = lon
                end
            end
            pooled = use_gram_isolated_pool && _gram_isolated_pool_batch_eval!(
                p.shared_buffers.densities,
                p.shared_buffers.temperatures,
                p.shared_buffers.winds,
                batch_model,
                alts,
                lats,
                lons,
                Float64(integrator.t),
                true,
                p;
                allotment_hint=decision.allotment
            )
            if !pooled
                getDensityBatch!(
                    p.shared_buffers.densities,
                    p.shared_buffers.temperatures,
                    p.shared_buffers.winds,
                    batch_model,
                    alts,
                    lats,
                    lons,
                    Float64(integrator.t),
                    true,
                    p
                )
            end
        elseif use_threads
            if isempty(density_models)
                density_model = fallback_density_model
                ParallelPolicy.threaded_foreach_persistent(:density_callback, num_sats, decision.allotment) do i
                    @inbounds update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
            else
                ParallelPolicy.threaded_foreach_persistent(:density_callback, num_sats, decision.allotment) do i
                    density_model = @inbounds density_models[i]
                    @inbounds update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
            end
        else
            if isempty(density_models)
                density_model = fallback_density_model
                @inbounds for i in 1:num_sats
                    update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
            else
                @inbounds for i in 1:num_sats
                    density_model = density_models[i]
                    update_density_sat!(i, p, u, integrator.t, segment_end_t, density_model, gram_density_caches)
                end
            end
        end
        if decision.policy_applied
            ParallelPolicy.record_policy_observation!(
                :density_callback;
                mode=decision.mode,
                num_items=num_sats,
                use_threads=use_threads,
                elapsed_ns=(time_ns() - started_ns)
            )
        end
    end

    return DiscreteCallback(condition, affect!, initialize=(cb, u, t, integrator) -> affect!(integrator))
end
