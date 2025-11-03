module Effectors

using LinearAlgebra, CSV, DataFrames
using ..Model
using ..Components
using ..Kinematics # We'll need this for rotate_to_body

export update_thrusters!, thrust_calculation_schmitt_trigger!, schmitt_trigger, integrate_impulse!

function update_thrusters!(link::Link, torque::AbstractVector{Float64}, t::Float64)
    """
    Updates the thrusters of the Link with the given torque vector.
    - `link`: The Link whose thrusters are to be updated.
    - `torque`: The torque vector to apply to the thrusters.
    """
    # Update the thruster Jacobian matrix to account for possible articulated joints
    link.J_thruster = MMatrix{3, length(link.thrusters), Float64}(zeros(3, length(link.thrusters))) # Reset the Jacobian matrix
    rot_to_body = rotate_to_body(link) # Get the rotation matrix to convert from inertial to body frame
    
    for (i, thruster) in enumerate(link.thrusters)
        normalize!(thruster.direction) # Ensure the thruster direction is a unit vector
        link.J_thruster[:, i] = cross(rot_to_body * thruster.location + link.r, rot_to_body * thruster.direction) # Update the Jacobian matrix with the r x F vector in the body frame
    end
    
    thrust_vector = pinv(link.J_thruster) * torque # Solve for the thrust vector using the Jacobian matrix
    thrust_vector .-= minimum(thrust_vector) # Ensure no negative thrust values
    
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
    ti = min(thrust / thruster.max_thrust * link.attitude_control_rate, link.attitude_control_rate) # Calculate the time interval for the thrust
    if ti < thruster.min_firing_time # If the time interval is less than the minimum firing time, set it to the minimum firing time or 0, based on schmitt trigger
        ti = schmitt_trigger(ti, thruster.level_on, thruster.level_off) * thruster.min_firing_time # Use Schmitt trigger to determine if the thruster should fire
    end
    
    CSV.write("thruster_debug.csv", DataFrame(time=time, on_time_request=ti, thrust_req=thrust), append=true)
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
    ω = thruster.cutoff_frequency # Get the cutoff frequency of the thruster
    κ = thruster.κ # Get the current thrust factor
    
    total_integrated_thrust = thruster.max_thrust * (on_time_request + (κ - 1) / ω * (1 - exp(-ω * on_time_request))) # Calculate the total impulse
    κ = 1 + (κ - 1) * exp(-ω * on_time_request) # Calculate the final thrust factor
    
    if on_time_request < link.attitude_control_rate
        total_integrated_thrust += thruster.max_thrust * κ / ω * (1 - exp(-ω * (link.attitude_control_rate - on_time_request))) # Add the impulse from ramp-down if applicable
        κ *= exp(-ω * (link.attitude_control_rate - on_time_request)) # Update the final thrust factor after ramp-down
    end
    
    thruster.κ = κ # Update the thrust factor in the thruster
    return total_integrated_thrust # Return the total impulse integrated over the control period
end

end # module Effectors