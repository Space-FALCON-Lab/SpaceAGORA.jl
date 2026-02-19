module SimulationCallbacks
include("../utils/Reference_system.jl") # Get the reference system types for the callback

using DifferentialEquations
using LinearAlgebra
using SPICE
using Dates
using ..EnvironmentModels: getDensity
using ..AbstractTypes: AbstractPlanet
export get_callbacks

function get_callbacks(num_sats::Int, effectors::Tuple)
    
    callbacks = CallbackSet(update_planet_frame_callback(), get_density_callback(num_sats), get_impact_callback(num_sats))
    if num_sats == 1
        callbacks = CallbackSet(callbacks, get_orbit_end_callback(), get_drag_state_callback())
    end
    return callbacks
end
function update_planet_frame_callback()
    condition(u, t, integrator) = true
    et_start = Ref{Float64}(NaN) # Store the start epoch in a Ref so it can be updated in the affect! function
    function affect!(integrator)
        p = integrator.p
        planet = p.args.environment_model.planet
        # initial_time = p.args.initial_time
        # current_epoch = from_utc(DateTime(initial_time.year, initial_time.month, initial_time.day, initial_time.hour, initial_time.minute, initial_time.second)) + integrator.t*seconds
        # utc_time = to_utc(current_epoch)
        # et = utc2et(utc_time)
        # planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(planet.name)", et))*planet.J2000_to_pci'
        # if isnan(et_start[])
        #     initial_time = p.args.initial_time
        #     start_epoch = from_utc(DateTime(
        #         initial_time.year,
        #         initial_time.month,
        #         initial_time.day,
        #         initial_time.hour,
        #         initial_time.minute,
        #         initial_time.second
        #     ))
        #     et_start[] = utc2et(to_utc(start_epoch))
        # end

        et = et_start[] + integrator.t
        planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(planet.name)", et)) * planet.J2000_to_pci'
    end

    function init_affect!(cb, u, t, integrator)
        p = integrator.p
        initial_time = p.args.initial_time
        start_epoch = from_utc(DateTime(
            initial_time.year,
            initial_time.month,
            initial_time.day,
            initial_time.hour,
            initial_time.minute,
            initial_time.second
        ))
        et_start[] = utc2et(to_utc(start_epoch))
        affect!(integrator) # Call the affect! function at the start of the simulation to initialize the planet frame
    end

    return DiscreteCallback(condition, affect!, initialize=init_affect!)
end
# Factory function to build the callback
function get_density_callback(num_sats::Int)
    
    # Logic for when the callback should trigger
    # Returning true means it runs at every integration step
    condition(u, t, integrator) = true

    # Logic for what happens when it triggers
    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        
        # Update shared buffers in the parameters
        @inbounds for i in 1:num_sats
            # Your density math here
            @views begin
                pos = SVector{3, Float64}(u.sc[i].pos) # Example of how to get the position of the satellite from the state vector (change to u.sc[i].r if using StructArrays in Complete_passage)
                rp, _ = r_intor_p!(pos, SVector{3, Float64}(u.sc[i].vel), p.args.environment_model.planet) # Convert position to planet-fixed frame for density lookup, also convert velocity if needed for density model
                alt,lat,lon = rtolatlong(rp, p.args.environment_model.planet) # Convert position to altitude, latitude, and longitude for density lookup
                ρ, T, wind_vec = getDensity(p.args.environment_model.density_model, alt, lat, lon, integrator.t, true, p) # Get the density from the density model and store in the shared buffer for this satellite
                p.shared_buffers.densities[i] = ρ
                p.shared_buffers.temperatures[i] = T
                p.shared_buffers.winds[i] = wind_vec
                # p.shared_buffers.densities[i], p.shared_buffers.temperatures[i], p.shared_buffers.winds[i] = getDensity(p.args.environment_model.density_model, alt, lat, lon, integrator.t, true, p) # Example of how to update the density buffer for each satellite
            end
        end
    end

    return DiscreteCallback(condition, affect!, initialize=(cb, u, t, integrator) -> affect!(integrator))
end

function get_impact_callback(num_sats::Int)
    function condition(u, t, integrator)
        p = integrator.p
        # Check for impacts based on the current state (placeholder logic, replace with actual impact condition)
        @inbounds for i in 1:num_sats
            if abs(norm(u.sc[i].pos) - p.args.environment_model.planet.Rp_e) < 50e3 # Example condition for impact (e.g., altitude below 80 km)
                return true
            end
        end
        return false
    end

    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        # Check for impacts and log them (placeholder logic)
        @inbounds for i in 1:num_sats
            if abs(norm(u.sc[i].pos) - p.args.environment_model.planet.Rp_e) < 50e3 # Example condition for impact (e.g., altitude below 80 km)
                println("Impact detected for satellite $i at time $(integrator.t) seconds!")
                p.is_active[i] = false # Mark the satellite as inactive after impact
                if all(p.is_active .== false) # If all satellites are inactive, we can stop the simulation
                    println("All satellites have impacted. Stopping simulation.")
                    terminate!(integrator)
                end
                # Log impact details to a file or data structure as needed
            end
        end
    end

    return DiscreteCallback(condition, affect!)
end

# Only call if simulating a single satellite
function get_orbit_end_callback()
    function condition(u, t, integrator)
        # Implement logic to determine if the satellite has completed an orbit (defined as reaching apoapsis again after passing periapsis)
        OE = rvtoorbitalelement(SVector{3, Float64}(u.sc[1].pos), SVector{3, Float64}(u.sc[1].vel), integrator.p.args.environment_model.planet)
        return OE[6] - pi
    end
    function affect!(integrator)
        p = integrator.p
        u = integrator.u
        p.orbit_counter += 1 # Increment the orbit counter in the shared buffers
        # Check for orbit completion and log it (placeholder logic)
        println("Orbit $(p.orbit_counter) completed at time $(integrator.t) seconds!")
    end

    return ContinuousCallback(condition, affect!, nothing)
end

# Only call if simulating single satellite
# Callback to adjust the max timestep depending on whether the satellite is in the atmosphere
function get_drag_state_callback()
    condition(u, t, integrator) = norm(u.sc[1].pos) - integrator.p.args.environment_model.planet.Rp_e - integrator.p.args.environment_model.EI*1e3
    function affect_upcrossing!(integrator)
        p = integrator.p
        u = integrator.u
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_orbit # Increase the maximum timestep when exiting the atmosphere
    end

    function affect_downcrossing!(integrator)
        p = integrator.p
        u = integrator.u
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_atmosphere # Decrease the maximum timestep when entering the atmosphere
    end

    return ContinuousCallback(condition, affect_upcrossing!, affect_downcrossing!)
end
end # module