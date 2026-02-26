module SimulationCallbacks
include("../utils/Reference_system.jl") # Get the reference system types for the callback

using DifferentialEquations
using LinearAlgebra
using SPICE
using Dates
using ..SimulationModel: SPICE_LOCK
using ..EnvironmentModels
using ..EnvironmentModels: getDensity, NoAtmosphereModel
using ..DynamicEffectors: BaseThrusterModel, AerodynamicCoefficientConstant, AerodynamicCoefficientfM, AerodynamicCoefficientNoBallisticFlight
using ..AbstractTypes: AbstractPlanet, AbstractDensityModel
using ..ConfigTypes: SaveData
using ..ControlEffectors: calcControlEffect!
using ..GuidanceEffectors: calcGuidanceEffect!
using ..SimConfig: SimulationConfiguration, MissionOrbits
export get_callbacks

@inline callback_verbose(integrator) = integrator.p.args.simulation_settings.verbose
@inline callback_use_invokelatest() = get(ENV, "SPACEAGORA_DEV_HOT_RELOAD", "0") == "1"

@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _density_callback_parallel_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_DENSITY_CALLBACK_PARALLEL", "auto")))
    if mode in ("off", "none", "serial", "0", "false", "no")
        return :off
    elseif mode in ("on", "thread", "threads", "1", "true", "yes")
        return :on
    elseif mode == "auto"
        return :auto
    end
    throw(ArgumentError("Unsupported SPACEAGORA_DENSITY_CALLBACK_PARALLEL='$mode'. Use one of: off, auto, on."))
end

@inline function _density_callback_thread_threshold()::Int
    raw = strip(get(ENV, "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD", "8"))
    threshold = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD must be an integer, got '$raw'"))
    end
    return max(1, threshold)
end

@inline function _control_callback_parallel_mode()::Symbol
    mode = lowercase(strip(get(ENV, "SPACEAGORA_CONTROL_CALLBACK_PARALLEL", "auto")))
    if mode in ("off", "none", "serial", "0", "false", "no")
        return :off
    elseif mode in ("on", "thread", "threads", "1", "true", "yes")
        return :on
    elseif mode == "auto"
        return :auto
    end
    throw(ArgumentError("Unsupported SPACEAGORA_CONTROL_CALLBACK_PARALLEL='$mode'. Use one of: off, auto, on."))
end

@inline function _control_callback_thread_threshold()::Int
    raw = strip(get(ENV, "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD", "8"))
    threshold = try
        parse(Int, raw)
    catch
        throw(ArgumentError("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD must be an integer, got '$raw'"))
    end
    return max(1, threshold)
end

# Extend this for custom user density models as needed:
# SimulationModel.SimulationCallbacks.density_model_threadsafe(::MyDensityModel) = true
@inline density_model_threadsafe(::AbstractDensityModel)::Bool = false
@inline density_model_threadsafe(::NoAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.ExponentialAtmosphereModel)::Bool = true
@inline density_model_threadsafe(::EnvironmentModels.PolynomialFitAtmosphereModel)::Bool = true

@inline function _density_callback_use_threads(args::SimulationConfiguration, num_sats::Int)::Bool
    if Threads.nthreads() <= 1 || num_sats <= 1
        return false
    end
    mode = _density_callback_parallel_mode()
    mode == :off && return false

    model = args.environment_model.density_model
    model_threadsafe = density_model_threadsafe(model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_ASSUME_THREADSAFE", false)
        return false
    end
    if mode == :on
        return true
    end
    return num_sats >= _density_callback_thread_threshold()
end

# Extend this for custom user control models as needed:
# SimulationModel.SimulationCallbacks.control_model_threadsafe(::MyControlModel) = true
@inline control_model_threadsafe(::Any)::Bool = false
@inline control_model_threadsafe(::BaseThrusterModel)::Bool = true

@inline function _control_callback_use_threads(control_model, num_sats::Int, use_invokelatest::Bool)::Bool
    if use_invokelatest || Threads.nthreads() <= 1 || num_sats <= 1
        return false
    end
    mode = _control_callback_parallel_mode()
    mode == :off && return false

    model_threadsafe = control_model_threadsafe(control_model)
    if !model_threadsafe && !_parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE", false)
        return false
    end
    if mode == :on
        return true
    end
    return num_sats >= _control_callback_thread_threshold()
end

@inline function _uses_atmospheric_dynamic_effector(effectors::Tuple)::Bool
    @inbounds for effector in effectors
        if effector isa AerodynamicCoefficientConstant || effector isa AerodynamicCoefficientfM || effector isa AerodynamicCoefficientNoBallisticFlight
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

@inline function _requires_drag_state_callback(effectors::Tuple, args::SimulationConfiguration)::Bool
    if !_requires_density_callback(effectors, args)
        return false
    end
    tol = args.integration_tolerances
    return tol.dt_max_atmosphere != tol.dt_max_orbit ||
           tol.reltol_atmosphere != tol.reltol_orbit ||
           tol.abstol_atmosphere != tol.abstol_orbit
end

function get_callbacks(num_sats::Int, effectors::Tuple, args::SimulationConfiguration)::CallbackSet
    callbacks = Any[
        get_impact_callback(num_sats),
        update_planet_frame_callback()
    ]

    if _requires_density_callback(effectors, args)
        push!(callbacks, get_density_callback(num_sats, args))
    end

    if _requires_orbit_end_callback(args)
        push!(callbacks, get_orbit_end_callback(num_sats))
    end

    if _requires_drag_state_callback(effectors, args)
        push!(callbacks, get_drag_state_callback(num_sats))
    end

    append!(callbacks, get_control_callbacks(num_sats, args))
    append!(callbacks, get_guidance_callbacks(num_sats, args))

    return CallbackSet(callbacks...)
end
function update_planet_frame_callback()
    condition(u, t, integrator) = true
    et_start = Ref{Float64}(0.0) # Store the start epoch in a Ref so it can be updated in the affect! function
    function affect!(integrator)
        p = integrator.p
        planet = p.args.environment_model.planet
        et = et_start[] + integrator.t
        lock(SPICE_LOCK) do
            planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_$(planet.name)", et)) * planet.J2000_to_pci'
        end
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
        lock(SPICE_LOCK) do
            et_start[] = utc2et(to_utc(start_epoch))
        end
        affect!(integrator) # Call the affect! function at the start of the simulation to initialize the planet frame
    end

    return DiscreteCallback(condition, affect!, initialize=init_affect!)
end
# Factory function to build the callback
function get_density_callback(num_sats::Int, args::SimulationConfiguration)
    use_threads = _density_callback_use_threads(args, num_sats)
    function update_density_sat!(i::Int, p, u, t)
        pos = SVector{3, Float64}(u.sc[i].pos)
        vel = SVector{3, Float64}(u.sc[i].vel)
        rp, _ = r_intor_p!(pos, vel, p.args.environment_model.planet)
        alt, lat, lon = rtolatlong(rp, p.args.environment_model.planet)
        ρ, T, wind_vec = getDensity(p.args.environment_model.density_model, alt, lat, lon, t, true, p)
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

        if use_threads
            Threads.@threads for i in 1:num_sats
                @inbounds update_density_sat!(i, p, u, integrator.t)
            end
        else
            @inbounds for i in 1:num_sats
                update_density_sat!(i, p, u, integrator.t)
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
                if callback_verbose(integrator)
                    println("Impact detected for satellite $i at time $(integrator.t) seconds!")
                end
                p.is_active[i] = false # Mark the satellite as inactive after impact
                if all(p.is_active .== false) # If all satellites are inactive, we can stop the simulation
                    if callback_verbose(integrator)
                        println("All satellites have impacted. Stopping simulation.")
                    end
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
        if callback_verbose(integrator)
            println("Orbit $(p.orbit_counter[idx]) completed by Satellite $idx at time $(integrator.t) seconds!")
        end
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
        if callback_verbose(integrator)
            println("Switching to space integration at time $(integrator.t) seconds!")
        end
        # sleep(3.0)
        integrator.opts.dtmax = p.args.integration_tolerances.dt_max_orbit # Increase the maximum timestep when exiting the atmosphere
        integrator.opts.reltol = p.args.integration_tolerances.reltol_orbit # Adjust the tolerances when exiting the atmosphere
        integrator.opts.abstol = p.args.integration_tolerances.abstol_orbit
    end

    function affect_downcrossing!(integrator, idx::Int64)
        p = integrator.p
        u = integrator.u
        if callback_verbose(integrator)
            println("Switching to atmosphere integration at time $(integrator.t) seconds!")
        end
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
    use_invokelatest = callback_use_invokelatest()
    callbacks = Vector{DiscreteCallback}(undef, length(guidance_models))
    for i in eachindex(guidance_models)
        guidance_model = guidance_models[i]
        guidance_rate = guidance_rates[i]
        # Implement a callback for this guidance model that triggers at the specified guidance rate and calculates the guidance commands based on the current state and the guidance model
        # The calculated guidance commands should be stored in the shared buffers for use in the dynamics calculations
        guidance_func = (integrator) -> begin
            if use_invokelatest
                # Dev mode: keep Revise/hot-reload workflows free of world-age errors.
                Base.invokelatest(calcGuidanceEffect!, guidance_model, integrator.u, integrator.p, integrator.t, i)
            else
                # Production mode: direct dispatch avoids invokelatest overhead.
                calcGuidanceEffect!(guidance_model, integrator.u, integrator.p, integrator.t, i)
            end
        end
        callbacks[i] = PeriodicCallback(guidance_func, guidance_rate)
    end
    return callbacks
end

function get_control_callbacks(num_sats::Int, args::SimulationConfiguration)::Vector{DiscreteCallback}
    # Perform the control effects' calculations at specific rates given by the control_rates field in the ControlModel
    control_models = args.control_model.control_effectors
    control_rates = args.control_model.control_rates
    use_invokelatest = callback_use_invokelatest()
    callbacks = Vector{DiscreteCallback}(undef, length(control_models))
    for i in eachindex(control_models)
        control_model = control_models[i]
        control_rate = control_rates[i]
        use_threads = _control_callback_use_threads(control_model, num_sats, use_invokelatest)
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
        apply_control! = if use_invokelatest
            # Dev mode: keep Revise/hot-reload workflows free of world-age errors.
            (integrator, sat_idx) -> Base.invokelatest(calcControlEffect!, control_model, integrator.u, integrator.p, integrator.t, sat_idx)
        else
            # Production mode: direct dispatch avoids invokelatest overhead.
            (integrator, sat_idx) -> calcControlEffect!(control_model, integrator.u, integrator.p, integrator.t, sat_idx)
        end
        control_func = (integrator) -> begin
            if use_threads
                Threads.@threads for sat_idx in 1:num_sats
                    @inbounds apply_control!(integrator, sat_idx)
                end
            else
                @inbounds for sat_idx in 1:num_sats
                    apply_control!(integrator, sat_idx)
                end
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

        if callback_verbose(integrator)
            println("Periapsis reached for Satellite $idx at time $(integrator.t) seconds!")
        end
        # Save the state at periapsis to a file or data structure as needed
    end

    return VectorContinuousCallback(condition!, affect!, nothing, num_sats)
end
end # module
