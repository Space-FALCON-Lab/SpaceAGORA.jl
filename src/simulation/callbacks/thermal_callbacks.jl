function get_thermal_callback(num_sats::Int, args::SimulationConfiguration)
    function update_thermal_sat!(i::Int, p, u)
        links = p.args.dynamics_model.spacecraft[i].links
        n_links = length(links)
        heat_rates = p.shared_buffers.heat_rates[i]
        if length(heat_rates) != n_links
            resize!(heat_rates, n_links)
        end
        fill!(heat_rates, 0.0)

        ρ = p.shared_buffers.densities[i]
        T = p.shared_buffers.temperatures[i]
        wind = p.shared_buffers.winds[i]
        if !isfinite(ρ) || !isfinite(T) || ρ <= 0.0 || T <= 0.0
            return nothing
        end

        planet = p.args.environment_model.planet
        thermal_model = p.args.environment_model.thermal_model
        pos_ii = SVector{3, Float64}(u.sc[i].pos)
        vel_ii = SVector{3, Float64}(u.sc[i].vel)
        pos_pp, vel_pp = r_intor_p!(pos_ii, vel_ii, planet)
        alt_lat_lon = rtolatlong(pos_pp, planet)
        uD, uN, uE = latlongtoNED(alt_lat_lon)
        wE, wN, wU = wind
        wind_pp = wN * uN + wE * uE - wU * uD
        vel_pp_rw = vel_pp + wind_pp
        v = norm(vel_pp_rw)
        sound_velocity = sqrt(planet.γ * planet.R * T)
        if !isfinite(v) || !isfinite(sound_velocity) || v <= 0.0 || sound_velocity <= 0.0
            return nothing
        end
        mach = v / sound_velocity
        S = sqrt(planet.γ * 0.5) * mach
        @inbounds for j in eachindex(links)
            α = links[j].α
            if !isfinite(α)
                continue
            end
            qdot = getHeatRate(thermal_model, S, T, ρ, v, α)
            heat_rates[j] = (isfinite(qdot) && qdot > 0.0) ? qdot : 0.0
        end
        return nothing
    end

    condition(u, t, integrator) = true

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        decision = _thermal_callback_thread_decision(num_sats)
        use_threads = decision.use_threads
        started_ns = time_ns()
        if use_threads
            ParallelPolicy.threaded_foreach_persistent(:thermal_callback, num_sats, decision.allotment) do i
                @inbounds update_thermal_sat!(i, p, u)
            end
        else
            @inbounds for i in 1:num_sats
                update_thermal_sat!(i, p, u)
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

