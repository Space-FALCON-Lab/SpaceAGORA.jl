module ThrusterHooks

using LinearAlgebra, CSV, DataFrames
using ..SpacecraftModels
using ..Components
using ..Kinematics

export thrust_calculation_schmitt_trigger!, schmitt_trigger, integrate_impulse!

# Thruster firing dynamics over a link's jets: the Schmitt trigger, impulse
# integration and the ramp factor. Torque allocation over mounted jets is
# `allocate_actuators` over `allocation_matrix` (the actuator layer).

@inline thruster_debug_enabled() = get(ENV, "SPACEAGORA_DEBUG_THRUSTER", "0") == "1"

function thrust_calculation_schmitt_trigger!(link::Link, thruster::Thruster, thrust::Float64, time::Float64;
                                             period::Float64=link.attitude_control_rate)
    """
    Applies a Schmitt trigger to the thrusters of the Link.
    """
    if !isfinite(thruster.max_thrust) || thruster.max_thrust <= 0.0
        thruster.thrust = 0.0
        return
    end
    ti = min(thrust / thruster.max_thrust * period, period) # Calculate the time interval for the thrust
    if ti < thruster.min_firing_time # If the time interval is less than the minimum firing time, set it to the minimum firing time or 0, based on schmitt trigger
        ti = schmitt_trigger(ti, thruster.level_on, thruster.level_off) * thruster.min_firing_time # Use Schmitt trigger to determine if the thruster should fire
    end
    
    if thruster_debug_enabled()
        CSV.write("thruster_debug.csv", DataFrame(time=time, on_time_request=ti, thrust_req=thrust), append=true)
    end
    total_integrated_thrust = integrate_impulse!(link, thruster, ti, time; period=period) # Integrate the impulse over the time interval
    
    thruster.thrust = total_integrated_thrust / period # Update the average thrust value in the thruster
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

function integrate_impulse!(link::Link, thruster::Thruster, on_time_request::Float64, time::Float64;
                            period::Float64=link.attitude_control_rate)
    """
    Integrates the impulse and thrust factor over the on-time request period.
    """
    ti = clamp(on_time_request, 0.0, period)
    ω = abs(thruster.cutoff_frequency) # Get the cutoff frequency of the thruster
    if !isfinite(ω) || ω < 1e-9
        ω = 1e-9
    end
    κ = thruster.κ # Get the current thrust factor
    
    one_minus_exp_on = -expm1(-ω * ti)
    total_integrated_thrust = thruster.max_thrust * (ti + (κ - 1) / ω * one_minus_exp_on) # Calculate the total impulse
    κ = 1 + (κ - 1) * exp(-ω * ti) # Calculate the final thrust factor
    
    if ti < period
        ramp_down_dt = period - ti
        one_minus_exp_down = -expm1(-ω * ramp_down_dt)
        total_integrated_thrust += thruster.max_thrust * κ / ω * one_minus_exp_down # Add the impulse from ramp-down if applicable
        κ *= exp(-ω * ramp_down_dt) # Update the final thrust factor after ramp-down
    end
    
    thruster.κ = κ # Update the thrust factor in the thruster
    return total_integrated_thrust # Return the total impulse integrated over the control period
end

end # module ThrusterHooks
