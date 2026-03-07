module ThrusterHooks

using LinearAlgebra, CSV, DataFrames
using ..SpacecraftModels
using ..Components
using ..Kinematics # We'll need this for rotate_to_body

export update_thrusters!, thrust_calculation_schmitt_trigger!, schmitt_trigger, integrate_impulse!

@inline thruster_debug_enabled() = get(ENV, "SPACEAGORA_DEBUG_THRUSTER", "0") == "1"

function update_thrusters!(link::Link, torque::AbstractVector{Float64}, t::Float64)
    """
    Updates the thrusters of the Link with the given torque vector.
    - `link`: The Link whose thrusters are to be updated.
    - `torque`: The torque vector to apply to the thrusters.
    """
    if isempty(link.thrusters)
        link.J_thruster = Matrix{Float64}(undef, 3, 0)
        return
    end

    # Update the thruster Jacobian matrix to account for possible articulated joints
    link.J_thruster = Matrix{Float64}(zeros(3, length(link.thrusters))) # Reset the Jacobian matrix
    rot_to_body = rotate_to_body(link) # Get the rotation matrix to convert from inertial to body frame
    
    for (i, thruster) in enumerate(link.thrusters)
        direction_norm = norm(thruster.direction)
        if !isfinite(direction_norm) || direction_norm == 0.0
            link.J_thruster[:, i] .= 0.0
            continue
        end
        thruster.direction ./= direction_norm # Ensure the thruster direction is a unit vector
        link.J_thruster[:, i] = cross(rot_to_body * thruster.location + link.r, rot_to_body * thruster.direction) # Update the Jacobian matrix with the r x F vector in the body frame
    end
    
    thrust_vector = pinv(link.J_thruster) * torque # Solve for the thrust vector using the Jacobian matrix
    if !all(isfinite, thrust_vector)
        thrust_vector .= 0.0
    else
        min_thrust = minimum(thrust_vector)
        if isfinite(min_thrust) && min_thrust < 0.0
            thrust_vector .-= min_thrust # Shift only when needed to avoid negative thrust values.
        end
        thrust_vector .= max.(thrust_vector, 0.0)
    end
    
    for (i, thruster) in enumerate(link.thrusters)
        thruster.thrust = thrust_vector[i] # Update the requested thrust in the thruster
        # Call the specific calculation function
        thrust_calculation_schmitt_trigger!(link, thruster, thrust_vector[i], t) 
    end
end

function thrust_calculation_schmitt_trigger!(link::Link, thruster::Thruster, thrust::Float64, time::Float64)
    """
    Applies a Schmitt trigger to the thrusters of the Link.
    """
    if !isfinite(thruster.max_thrust) || thruster.max_thrust <= 0.0
        thruster.thrust = 0.0
        return
    end
    ti = min(thrust / thruster.max_thrust * link.attitude_control_rate, link.attitude_control_rate) # Calculate the time interval for the thrust
    if ti < thruster.min_firing_time # If the time interval is less than the minimum firing time, set it to the minimum firing time or 0, based on schmitt trigger
        ti = schmitt_trigger(ti, thruster.level_on, thruster.level_off) * thruster.min_firing_time # Use Schmitt trigger to determine if the thruster should fire
    end
    
    if thruster_debug_enabled()
        CSV.write("thruster_debug.csv", DataFrame(time=time, on_time_request=ti, thrust_req=thrust), append=true)
    end
    total_integrated_thrust = integrate_impulse!(link, thruster, ti, time) # Integrate the impulse over the time interval
    
    thruster.thrust = total_integrated_thrust / link.attitude_control_rate # Update the average thrust value in the thruster
end

function schmitt_trigger(input::Float64, level_on::Float64, level_off::Float64)
    """
    Implements a Schmitt trigger.
    """
    state = 0.0
    if input > level_on
        state = 1.0
    elseif input < level_off
        state = 0.0
    end
    return state
end

function integrate_impulse!(link::Link, thruster::Thruster, on_time_request::Float64, time::Float64)
    """
    Integrates the impulse and thrust factor over the on-time request period.
    """
    ti = clamp(on_time_request, 0.0, link.attitude_control_rate)
    ω = abs(thruster.cutoff_frequency) # Get the cutoff frequency of the thruster
    if !isfinite(ω) || ω < 1e-9
        ω = 1e-9
    end
    κ = thruster.κ # Get the current thrust factor
    
    one_minus_exp_on = -expm1(-ω * ti)
    total_integrated_thrust = thruster.max_thrust * (ti + (κ - 1) / ω * one_minus_exp_on) # Calculate the total impulse
    κ = 1 + (κ - 1) * exp(-ω * ti) # Calculate the final thrust factor
    
    if ti < link.attitude_control_rate
        ramp_down_dt = link.attitude_control_rate - ti
        one_minus_exp_down = -expm1(-ω * ramp_down_dt)
        total_integrated_thrust += thruster.max_thrust * κ / ω * one_minus_exp_down # Add the impulse from ramp-down if applicable
        κ *= exp(-ω * ramp_down_dt) # Update the final thrust factor after ramp-down
    end
    
    thruster.κ = κ # Update the thrust factor in the thruster
    return total_integrated_thrust # Return the total impulse integrated over the control period
end

end # module ThrusterHooks
