function get_impact_callback(num_sats::Int)
    function condition!(out, u, t, integrator)
        p = integrator.p
        Rp_e = p.args.environment_model.planet.Rp_e
        @inbounds for i in 1:num_sats
            out[i] = norm(u.sc[i].pos) - Rp_e
        end
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        if p.is_active[idx]
            if callback_verbose(integrator)
                println("Impact detected for satellite $idx at time $(integrator.t) seconds!")
            end
            p.is_active[idx] = false
            if all(p.is_active .== false)
                if callback_verbose(integrator)
                    println("All satellites have impacted. Stopping simulation.")
                end
                terminate!(integrator)
            end
        end
    end

    return VectorContinuousCallback(condition!, nothing, affect_downcrossing!, num_sats)
end

# Only call if simulating a single satellite
function get_orbit_end_callback(num_sats::Int)
    function condition!(out, u, t, integrator)
        # Use radial-velocity root events for orbit bookkeeping.
        # At apoapsis, dot(r,v) crosses + -> -, so -dot(r,v) crosses - -> + (upcrossing),
        # which matches the single affect! handler below.
        @inbounds for i in 1:num_sats
            pos = SVector{3, Float64}(u.sc[i].pos)
            vel = SVector{3, Float64}(u.sc[i].vel)
            out[i] = -dot(pos, vel)
        end
    end

    function affect!(integrator, idx::Int64)
        p = integrator.p

        p.orbit_counter[idx] += 1 # Increment the orbit counter in the shared buffers
        completed_orbits = p.orbit_counter[idx] - 1
        target_orbits = p.args.mission_configuration.number_of_orbits
        # Check for orbit completion and log it (placeholder logic)
        if callback_verbose(integrator)
            println("Orbit $(completed_orbits) completed by Satellite $idx at time $(integrator.t) seconds!")
        end

        # MissionOrbits should end when every active satellite reaches the requested
        # number of completed orbits, analogous to drag-passage event termination.
        if p.args.mission_configuration.mission_type == MissionOrbits && completed_orbits >= target_orbits
            all_active_reached_target = true
            @inbounds for sat_idx in eachindex(p.orbit_counter)
                if p.is_active[sat_idx] && (p.orbit_counter[sat_idx] - 1) < target_orbits
                    all_active_reached_target = false
                    break
                end
            end
            if all_active_reached_target
                if callback_verbose(integrator)
                    println("Target orbit count reached for all active satellites. Stopping simulation.")
                end
                if applicable(terminate!, integrator)
                    terminate!(integrator)
                end
            end
        end
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
end

function get_entry_end_callback(num_sats::Int, args::SimulationConfiguration)
    target_entries = _entry_target_count()
    target_entries > 0 || throw(ArgumentError("Entry target callback requires SPACEAGORA_ENTRY_TARGET_COUNT > 0"))
    entry_interface_m = args.environment_model.EI * 1e3
    entry_counter = zeros(Int64, num_sats)

    function condition!(out, u, t, integrator)
        p = integrator.p
        planet = p.args.environment_model.planet
        @inbounds for i in 1:num_sats
            if !p.is_active[i]
                out[i] = 1.0
                continue
            end
            alt = norm(u.sc[i].pos) - planet.Rp_e
            out[i] = alt - entry_interface_m
        end
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        if !p.is_active[idx]
            return nothing
        end

        entry_counter[idx] += 1
        completed_entries = entry_counter[idx]
        if callback_verbose(integrator)
            println("Entry $(completed_entries) detected for Satellite $idx at time $(integrator.t) seconds!")
        end

        if completed_entries >= target_entries
            all_active_reached_target = true
            @inbounds for sat_idx in eachindex(entry_counter)
                if p.is_active[sat_idx] && entry_counter[sat_idx] < target_entries
                    all_active_reached_target = false
                    break
                end
            end
            if all_active_reached_target
                if callback_verbose(integrator)
                    println("Target entry count reached for all active satellites. Stopping simulation.")
                end
                if applicable(terminate!, integrator)
                    terminate!(integrator)
                end
            end
        end
        return nothing
    end

    return VectorContinuousCallback(condition!, nothing, affect_downcrossing!, num_sats)
end

# Only call if simulating single satellite
# Callback to adjust the max timestep depending on whether the satellite is in the atmosphere
function get_drag_state_callback(num_sats::Int)
    # out = zeros(num_sats) # Output array to store the condition for each satellite
    condition!(out, u, t, integrator) = begin
        @inbounds for i in 1:length(u.sc)
            alt = norm(u.sc[i].pos) - integrator.p.args.environment_model.planet.Rp_e
            # println("Satellite $i altitude: $(alt) meters at time $(integrator.t) seconds")
            out[i] = alt - integrator.p.args.environment_model.EI*1e3 # Positive when above the atmosphere, negative when in the atmosphere
        end
        # return out
    end
    function affect_upcrossing!(integrator, idx::Int64)
        p = integrator.p
        if callback_verbose(integrator)
            println("Switching to space integration at time $(integrator.t) seconds!")
        end
        # sleep(3.0)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_orbit # Increase the maximum timestep when exiting the atmosphere
        reltol_new, abstol_new = _callback_tolerances_for_phase(
            integrator.opts.reltol,
            integrator.opts.abstol,
            p.args,
            false
        )
        integrator.opts.reltol = reltol_new # Adjust tolerances when exiting the atmosphere
        integrator.opts.abstol = abstol_new
        schedule_event_driven_thruster_controls!(integrator, idx)
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        if callback_verbose(integrator)
            println("Switching to atmosphere integration at time $(integrator.t) seconds!")
        end
        # sleep(3.0)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_atmosphere # Decrease the maximum timestep when entering the atmosphere
        reltol_new, abstol_new = _callback_tolerances_for_phase(
            integrator.opts.reltol,
            integrator.opts.abstol,
            p.args,
            true
        )
        integrator.opts.reltol = reltol_new # Adjust tolerances when entering the atmosphere
        integrator.opts.abstol = abstol_new
    end

    return VectorContinuousCallback(condition!, affect_upcrossing!, affect_downcrossing!, num_sats)
end

function get_quaternion_projection_callback(num_sats::Int, args::SimulationConfiguration)
    tol_q = max(1e-12, args.integration_tolerances.abstol_quaternion)
    condition(u, t, integrator) = begin
        p = integrator.p
        @inbounds for i in 1:num_sats
            if !p.is_active[i]
                continue
            end
            q = u.sc[i].q
            qnorm2 = dot(q, q)
            if !isfinite(qnorm2) || abs(qnorm2 - 1.0) > tol_q
                return true
            end
        end
        return false
    end

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        corrected = false
        @inbounds for i in 1:num_sats
            if !p.is_active[i]
                continue
            end
            q = u.sc[i].q
            qnorm2 = dot(q, q)
            if !(isfinite(qnorm2) && qnorm2 > eps(Float64))
                # Keep attitude state finite if a pathological quaternion is encountered.
                q .= (0.0, 0.0, 0.0, 1.0)
                corrected = true
                continue
            end
            if abs(qnorm2 - 1.0) > tol_q
                q ./= sqrt(qnorm2)
                corrected = true
            end
        end
        if corrected && callback_verbose(integrator)
            println("Quaternion projection applied at time $(integrator.t) seconds.")
        end
    end

    return DiscreteCallback(
        condition,
        affect!,
        initialize=(cb, u, t, integrator) -> affect!(integrator)
    )
end

function get_data_saving_callback(
    num_sats::Int,
    args::SimulationConfiguration,
    save_fields,
    saved_values=nothing
)
    saved_values = isnothing(saved_values) ? SavedValues(Float64, SaveData) : saved_values
    function save_func(u, t, integrator)
        return _save_snapshot(save_fields, u, t, integrator)
    end
    data_rate = args.mission_configuration.data_rate
    data_rate > 0.0 || throw(ArgumentError("mission_configuration.data_rate must be > 0.0, got $data_rate."))
    return SavingCallback(save_func, saved_values; saveat=data_rate, save_everystep=false)
end

@inline function _progress_interval_s()::Float64
    return max(0.0, _parse_float_env("SPACEAGORA_PROGRESS_INTERVAL_S", 0.0))
end

function get_progress_callback(num_sats::Int, args::SimulationConfiguration)
    interval_s = _progress_interval_s()
    interval_s > 0.0 || throw(ArgumentError("Progress callback requires SPACEAGORA_PROGRESS_INTERVAL_S > 0."))
    mission_end_s = args.mission_configuration.mission_time

    function progress_affect!(integrator)
        p = integrator.p
        t_s = Float64(integrator.t)
        percent = clamp(100.0 * t_s / mission_end_s, 0.0, 100.0)
        min_altitude_m = Inf
        lead_speed_m_s = 0.0

        @inbounds for sat_idx in 1:num_sats
            if !p.is_active[sat_idx]
                continue
            end
            pos = Float64.(integrator.u.sc[sat_idx].pos)
            vel = Float64.(integrator.u.sc[sat_idx].vel)
            alt_m = norm(pos) - p.args.environment_model.planet.Rp_e
            min_altitude_m = min(min_altitude_m, alt_m)
            if sat_idx == 1
                lead_speed_m_s = norm(vel)
            end
        end

        min_altitude_km = isfinite(min_altitude_m) ? min_altitude_m / 1e3 : NaN
        lead_speed_km_s = lead_speed_m_s / 1e3
        println(
            "Progress: t=$(round(t_s / 3600.0, digits=3)) hr / $(round(mission_end_s / 3600.0, digits=3)) hr " *
            "($(round(percent, digits=1))%), min_altitude=$(round(min_altitude_km, digits=3)) km, " *
            "speed=$(round(lead_speed_km_s, digits=3)) km/s"
        )
        return nothing
    end

    return PeriodicCallback(progress_affect!, interval_s)
end


function get_periapsis_save_callback(num_sats::Int)
    # Implement a callback to save the state of the simulation at periapsis for each orbit, which can be useful for analyzing the changes in the orbit after each pass through the atmosphere
    function condition!(out, u, t, integrator)
        @inbounds for i in 1:num_sats
            planet = integrator.p.args.environment_model.planet
            OE = rvtoorbitalelement(SVector{3, Float64}(u.sc[i].pos), SVector{3, Float64}(u.sc[i].vel), planet)
            out[i] = OE[6] # Return the true anomaly (ν) which is zero at periapsis
        end
    end

    function affect!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u

        if callback_verbose(integrator)
            println("Periapsis reached for Satellite $idx at time $(integrator.t) seconds!")
        end
        # Save the state at periapsis to a file or data structure as needed
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
end
