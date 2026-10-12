@inline function _heat_rate_buffer_for_sat!(p, sat_idx::Int)
    links = p.args.dynamics_model.spacecraft[sat_idx].links
    n_links = length(links)
    heat_rates = p.shared_buffers.heat_rates[sat_idx]
    if length(heat_rates) != n_links
        resize!(heat_rates, n_links)
    end
    fill!(heat_rates, 0.0)
    return heat_rates
end

function _compute_stage_heat_rates!(
    p,
    x,
    sat_idx::Int,
    t::Float64;
    use_buffered_density::Bool=false,
)
    links = p.args.dynamics_model.spacecraft[sat_idx].links
    isempty(links) && return _heat_rate_buffer_for_sat!(p, sat_idx)
    heat_rates = _heat_rate_buffer_for_sat!(p, sat_idx)
    engine = _simulation_engine_module()
    planet_frame = engine.sample_planet_frame(x, p, sat_idx, t)
    atmosphere = use_buffered_density ?
        engine.sample_buffered_atmosphere(x, p, sat_idx, t) :
        engine.sample_atmosphere(x, p, sat_idx, t; write_buffers=false)

    rho = atmosphere.rho_kg_m3
    T = atmosphere.temperature_k
    wind = atmosphere.wind_pp
    if !isfinite(rho) || !isfinite(T) || rho <= 0.0 || T <= 0.0
        return heat_rates
    end

    planet = p.args.environment_model.planet
    thermal_model = p.args.environment_model.thermal_model
    uD, uN, uE = latlongtoNED((planet_frame.alt_m, planet_frame.lat_rad, planet_frame.lon_rad))
    wE, wN, wU = wind
    wind_pp = wN * uN + wE * uE - wU * uD
    vel_pp_rw = planet_frame.vel_pp - wind_pp
    v = norm(vel_pp_rw)
    sound_velocity = sqrt(planet.γ * planet.R * T)
    if !isfinite(v) || !isfinite(sound_velocity) || v <= 0.0 || sound_velocity <= 0.0
        return heat_rates
    end

    mach = v / sound_velocity
    S = sqrt(planet.γ * 0.5) * mach
    @inbounds for j in eachindex(links)
        alpha = links[j].α
        if !isfinite(alpha)
            continue
        end
        qdot = getHeatRate(thermal_model, S, T, rho, v, alpha)
        heat_rates[j] = (isfinite(qdot) && qdot > 0.0) ? qdot : 0.0
    end
    return heat_rates
end

function get_thermal_callback(num_sats::Int, args::SimulationConfiguration)
    function update_thermal_sat!(i::Int, p, u, t::Float64)
        _compute_stage_heat_rates!(p, u.sc[i], i, t; use_buffered_density=true)
        return nothing
    end

    condition(u, t, integrator) = true

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        decision = _thermal_callback_thread_decision(p, num_sats)
        use_threads = decision.use_threads
        started_ns = time_ns()
        if use_threads
            ParallelPolicy.threaded_foreach_persistent(:thermal_callback, num_sats, decision.allotment) do i
                @inbounds update_thermal_sat!(i, p, u, Float64(integrator.t))
            end
        else
            @inbounds for i in 1:num_sats
                update_thermal_sat!(i, p, u, Float64(integrator.t))
            end
        end
        if decision.policy_applied
            ParallelPolicy.record_policy_observation!(
                :thermal_callback;
                mode=decision.mode,
                num_items=num_sats,
                use_threads=use_threads,
                elapsed_ns=(time_ns() - started_ns)
            )
        end
    end

    return DiscreteCallback(condition, affect!, initialize=(cb, u, t, integrator) -> affect!(integrator))
end
