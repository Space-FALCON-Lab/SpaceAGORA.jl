module SimulationCallbacks
include("../utils/Reference_system.jl") # Get the reference system types for the callback

using DifferentialEquations
using LinearAlgebra
using SPICE
using Dates
using ..EnvironmentModels: getDensity
using ..DynamicEffectors: BaseThrusterModel
using ..AbstractTypes: AbstractPlanet
using ..ConfigTypes: SaveData
using ..ControlEffectors: calcControlEffect!
using ..GuidanceEffectors: calcGuidanceEffect!
using ..SimConfig: SimulationConfiguration
export get_callbacks

function get_callbacks(num_sats::Int, effectors::Tuple, args::SimulationConfiguration)::CallbackSet
    
    callbacks = CallbackSet(get_impact_callback(num_sats),
                            update_planet_frame_callback(), 
                            get_density_callback(num_sats), 
                            get_orbit_end_callback(num_sats),
                            get_drag_state_callback(num_sats),
                            get_control_callbacks(num_sats, args)...,
                            get_guidance_callbacks(num_sats, args)...)
     
    return callbacks
end
function update_planet_frame_callback()
    condition(u, t, integrator) = true
    et_start = Ref{Float64}(0.0) # Store the start epoch in a Ref so it can be updated in the affect! function
    function affect!(integrator)
        p = integrator.p
        planet = p.args.environment_model.planet
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
        # Check for impacts based on the current state
        @inbounds for i in 1:num_sats
            if p.is_active[i] && abs(norm(u.sc[i].pos) - p.args.environment_model.planet.Rp_e) < 50e3 # Example condition for impact (e.g., altitude below 80 km)
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
            if p.is_active[i] && abs(norm(u.sc[i].pos) - p.args.environment_model.planet.Rp_e) < 50e3 # Example condition for impact (e.g., altitude below 80 km)
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
function get_orbit_end_callback(num_sats::Int)
    function condition!(out, u, t, integrator)
        # if integrator.sol.retcode != :Default
        #     return false
        # end
        # Implement logic to determine if the satellite has completed an orbit (defined as reaching apoapsis again after passing periapsis)
        @inbounds for i in 1:num_sats
            OE = rvtoorbitalelement(SVector{3, Float64}(u.sc[i].pos), SVector{3, Float64}(u.sc[i].vel), integrator.p.args.environment_model.planet)
            out[i] = OE[6] - pi # Return the angle of periapsis minus pi (i.e., when the satellite is at apoapsis)
        end
    end

    function affect!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u

        p.orbit_counter[idx] += 1 # Increment the orbit counter in the shared buffers
        # Check for orbit completion and log it (placeholder logic)
        println("Orbit $(p.orbit_counter[idx]) completed by Satellite $idx at time $(integrator.t) seconds!")
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
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
        u = integrator.u
        println("Switching to space integration at time $(integrator.t) seconds!")
        # sleep(3.0)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_orbit # Increase the maximum timestep when exiting the atmosphere
        integrator.opts.reltol = p.args.integration_tolerances.reltol_orbit # Adjust the tolerances when exiting the atmosphere
        integrator.opts.abstol = p.args.integration_tolerances.abstol_orbit
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u
        println("Switching to atmosphere integration at time $(integrator.t) seconds!")
        # sleep(3.0)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_atmosphere # Decrease the maximum timestep when entering the atmosphere
        integrator.opts.reltol = p.args.integration_tolerances.reltol_atmosphere # Adjust the tolerances when entering the atmosphere
        integrator.opts.abstol = p.args.integration_tolerances.abstol_atmosphere
    end

    return VectorContinuousCallback(condition!, affect_upcrossing!, affect_downcrossing!, num_sats)
end

function get_data_saving_callback(num_sats::Int)
    # Implement a callback to save the state of the simulation at regular intervals defined by num_steps_to_save in the mission configuration
    # This will involve storing the state in the shared buffers and writing to a file when the buffer is full, then clearing the buffer for the next set of data
    # The condition can be based on the number of steps taken or the simulation time, depending on how you want to structure it
    saved_values = SaveValues(Float64, SaveData)
    function save_func(u, t, integrator)
        p = integrator.p

        # Save the current state to the SaveData struct

    end
    return SavingCallback(save_func, saved_values) # Placeholder, implement the actual logic for the data saving callback
end

function get_guidance_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    # Implement a callback to calculate guidance commands at each time step based on the current state and the guidance model defined in the simulation configuration
    guidance_models = args.guidance_model.guidance_effectors
    guidance_rates = args.guidance_model.guidance_rates
    callbacks = Vector{DiscreteCallback}(undef, length(guidance_models))
    for i in eachindex(guidance_models)
        guidance_model = guidance_models[i]
        guidance_rate = guidance_rates[i]
        # Implement a callback for this guidance model that triggers at the specified guidance rate and calculates the guidance commands based on the current state and the guidance model
        # The calculated guidance commands should be stored in the shared buffers for use in the dynamics calculations
        guidance_func = (integrator) -> begin
            # if integrator.sol.retcode == :Default
                calcGuidanceEffect!(guidance_model, integrator.u, integrator.p, integrator.t, i)
            # end
        end
        callbacks[i] = PeriodicCallback(guidance_func, guidance_rate)
    end
    return callbacks
end

function get_control_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    # Perform the control effects' calculations at specific rates given by the control_rates field in the ControlModel
    control_models = args.control_model.control_effectors
    control_rates = args.control_model.control_rates
    callbacks = Vector{DiscreteCallback}(undef, length(control_models))
    for i in eachindex(control_models)
        control_model = control_models[i]
        control_rate = control_rates[i]
        if control_model isa BaseThrusterModel
            n_slots = length(control_model.thrust)
            if n_slots != num_sats
                throw(ArgumentError(
                    "BaseThrusterModel vector length ($n_slots) must match number of spacecraft ($num_sats). " *
                    "Use one shared model with per-spacecraft vectors."
                ))
            end
        end
        # Each control effector callback runs at its own rate and updates
        # all spacecraft states. The spacecraft index is passed explicitly
        # to avoid conflating effector-index with spacecraft-index.
        control_func = (integrator) -> begin
            @inbounds for sat_idx in 1:num_sats
                calcControlEffect!(control_model, integrator.u, integrator.p, integrator.t, sat_idx)
            end
        end
        callbacks[i] = PeriodicCallback(control_func, control_rate)
    end
    return callbacks
end

function get_periapsis_save_callback(num_sats::Int)
    # Implement a callback to save the state of the simulation at periapsis for each orbit, which can be useful for analyzing the changes in the orbit after each pass through the atmosphere
    function condition!(out, u, t, integrator)
        @inbounds for i in 1:num_sats
            OE = rvtoorbitalelement(SVector{3, Float64}(u.sc[i].pos), SVector{3, Float64}(u.sc[i].vel), integrator.p.args.environment_model.planet)
            out[i] = OE[6] # Return the true anomaly (ν) which is zero at periapsis
        end
    end

    function affect!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u

        println("Periapsis reached for Satellite $idx at time $(integrator.t) seconds!")
        # Save the state at periapsis to a file or data structure as needed
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
end
end # module
