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
