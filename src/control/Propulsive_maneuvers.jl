include("../utils/Reference_system.jl")
using ComponentArrays

"""
calcControlForceTorque(controlModel::BaseThrusterModel, x::AbstractVector{Float64}, p::ODEParams, i::Int64, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}

Calculate the control force and torque based on the thruster model and current state, called in the dynamics loop to get the current thruster force
- `controlModel`: The thruster model containing thrust magnitudes, directions, burn times, and specific impulses for each thruster
- `x`: The current state vector of the spacecraft
- `p`: The ODE parameters containing simulation configuration and other relevant data
- `i`: The index of the spacecraft for which to calculate the control force/torque
- `t`: The current time in the simulation

Returns a tuple containing the total control force and torque as 3D vectors
"""
function calcControlForceTorque(controlModel::BaseThrusterModel, u::AbstractVector{Float64}, p::ODEParams, i::Int64, t::Float64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    # Calculate the control force and torque based on the thruster model and current state
    if i < 1 || i > length(controlModel.start_burn_time)
        return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    start_time = controlModel.start_burn_time[i]
    stop_time = controlModel.stop_burn_time[i]
    if t >= start_time && t <= stop_time
        thrust_mag = controlModel.thrust[i]
        vel_vec = SVector{3, Float64}(u.vel)
        vel_mag = norm(vel_vec)
        if vel_mag == 0.0 || !isfinite(vel_mag)
            return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
        end
        thrust_dir = normalize(vel_vec) # Prograde direction
        thrust_dir *= cos(controlModel.direction[i]) >= 0.0 ? 1.0 : -1.0 # 0 -> prograde, π -> retrograde
        force = thrust_mag * thrust_dir
        # For simplicity, assume torque is zero for now (can be updated later to include offset thrusters, gimbaled thrusters, etc.)
        torque = SVector{3, Float64}(0.0, 0.0, 0.0)
    else
        force = SVector{3, Float64}(0.0, 0.0, 0.0)
        torque = SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    return force, torque
end

"""
calcControlEffect!(controlModel::BaseThrusterModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
Calculate the control effect (force and torque) based on the control model and current state, and store it in the shared buffers for use in the dynamics calculations

Args
- `controlModel`: The thruster model containing thrust magnitudes, directions, burn times, and specific impulses for each thruster
- `u`: The current state vector of the spacecraft as a ComponentVector
- `p`: The ODE parameters containing simulation configuration and other relevant data
- `t`: The current time in the simulation
- `i`: The index of the spacecraft for which to calculate the control effect

Returns
- Updates the control force and torque in the shared buffers for the specified spacecraft index
"""
function calcControlEffect!(controlModel::BaseThrusterModel, u::ComponentVector, p::ODEParams, t::Float64, i::Int64)
    # Calculate the control effect (force and torque) based on the control model and current state, and store it in the shared buffers for use in the dynamics calculations
    if i < 1 || i > length(controlModel.start_burn_time)
        return
    end
    # Default behavior: track/recenter the burn window until ignition,
    # then lock it once the burn has started.
    start_time = controlModel.start_burn_time[i]
    stop_time = controlModel.stop_burn_time[i]
    if isfinite(start_time) && isfinite(stop_time) && stop_time > start_time && t >= start_time - 1e-9
        return
    end
    pos = SVector{3, Float64}(u.sc[i].pos)
    vel = SVector{3, Float64}(u.sc[i].vel)
    mass = u.sc[i].mass
    if !isfinite(mass) || mass <= 0.0
        return
    end
    # Calculate the current orbital elements from the state vector
    oe = try
        rvtoorbitalelement(pos, vel, p.args.environment_model.planet)
    catch
        return
    end
    a, e, _, _, _, ν = oe
    if !isfinite(a) || !isfinite(e) || !isfinite(ν) || a <= 0.0 || e < 0.0 || e >= 1.0
        return
    end
    
    # Don't start maneuver until spacecraft has exited the atmosphere and is past periapsis.
    alt = norm(pos) - p.args.environment_model.planet.Rp_e
    if alt >= p.args.environment_model.EI * 1000 - 1e-6 && ν < π - 1e-12
        # Calculate the burn time required to achieve the desired Δv based on the current mass and thrust
        Δv = controlModel.Δv[i]
        thrust_mag = controlModel.thrust[i]
        # Calculate the burn duration for a constant-thrust impulsive approximation.
        if thrust_mag <= 0.0 || !isfinite(thrust_mag)
            return
        end
        burn_time = mass * Δv / thrust_mag
        if !isfinite(burn_time) || burn_time < 0.0
            return
        end
        # println("e: $e, a: $a, Δv: $Δv, burn_time: $burn_time")
        # Estimate the time of apoapsis and set the burn to be symmetric about that time
        ψ = 2*atan(sqrt((1-e)/(1+e))*tan(ν/2))

        M = ψ - e*sin(ψ)
        n = sqrt(p.args.environment_model.planet.μ/a^3)
        if !isfinite(M) || !isfinite(n) || n <= 0.0
            return
        end
        time_to_apoapsis = (π - M)/n
        apoapsis_time = t + time_to_apoapsis
        if !isfinite(apoapsis_time)
            return
        end
        # Calculate start/end time as symmetric about the apoapsis time
        start_burn_time = apoapsis_time - burn_time/2
        stop_burn_time = apoapsis_time + burn_time/2
        
        # Update the start/end time fields in the control model for the current spacecraft
        controlModel.start_burn_time[i] = start_burn_time
        controlModel.stop_burn_time[i] = stop_burn_time
    end
end
