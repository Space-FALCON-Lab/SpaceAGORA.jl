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

@inline function _planet_lpi_at(p, t::Float64)::SMatrix{3, 3, Float64}
    planet = p.args.environment_model.planet
    ephemerides_model = p.args.environment_model.ephemerides_model
    et = p.shared_buffers.et_start[] + t
    pxform_counter = p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls
    cache_entry = p.shared_buffers.planet_frame_ephemeris_cache[]
    return if cache_entry isa PlanetFrameEphemerisCache
        cached = _planet_lpi_from_cache(cache_entry, et)
        cached === nothing ? _planet_lpi_from_backend(planet, ephemerides_model, et, pxform_counter) : cached
    else
        _planet_lpi_from_backend(planet, ephemerides_model, et, pxform_counter)
    end
end

@inline function _planet_relative_state(
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    planet,
    l_pi::SMatrix{3, 3, Float64}
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_pp = SVector{3, Float64}(l_pi * pos_ii)
    # planet.ω is the spin vector in the planet-fixed frame (pole = +z there), so
    # the transport term must be evaluated AFTER rotating into that frame; taking
    # the cross product in J2000 axes points the co-rotation velocity along the
    # J2000 pole instead of the body pole (~250 m/s misdirected at Mars periapsis).
    vel_pp = SVector{3, Float64}(l_pi * vel_ii - cross(planet.ω, pos_pp))
    return pos_pp, vel_pp
end

function update_planet_frame_callback()
    condition(u, t, integrator) = true
    function affect!(integrator)
        p = integrator.p
        p.args.environment_model.planet.L_PI .= _planet_lpi_at(p, integrator.t)
        # Invalidate _rhs_execution_plan's per-step cache (setup.jl / runtime_types.jl's
        # rhs_plan_step_cache) here since this callback already runs exactly once per
        # accepted step, unconditionally, for every non-backbone-mode solve -- the same
        # boundary the cache needs, with no extra CallbackSet entry required. No-op
        # (single false check) when SPACEAGORA_RHS_PLAN_STEP_CACHE is off (default).
        if _simulation_engine_module()._rhs_plan_step_cache_enabled()
            p.shared_buffers.rhs_plan_step_cache[] = nothing
        end
    end

    function init_affect!(cb, u, t, integrator)
        p = integrator.p
        p.shared_buffers.et_start[] = ephemerides_time_seconds(p.args.initial_time, p.args.environment_model.ephemerides_model)
        affect!(integrator) # Call the affect! function at the start of the simulation to initialize the planet frame
    end

    return DiscreteCallback(condition, affect!, initialize=init_affect!)
end
