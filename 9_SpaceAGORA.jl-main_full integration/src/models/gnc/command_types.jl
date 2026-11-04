module CommandTypes

export PropulsiveManeuverCommand, PropulsiveBurnPlan, AerobrakingControlCommand

Base.@kwdef struct PropulsiveManeuverCommand
    valid::Bool = false
    delta_v_mps::Float64 = 0.0
    direction_rad::Float64 = 0.0
    source_orbit::Int64 = 0
end

Base.@kwdef struct PropulsiveBurnPlan
    valid::Bool = false
    delta_v_mps::Float64 = 0.0
    direction_rad::Float64 = 0.0
    source_orbit::Int64 = 0
    start_burn_s::Float64 = -1.0
    stop_burn_s::Float64 = -1.0
    thrust_n::Float64 = 0.0
    isp_s::Float64 = 0.0
    commanded_impulse_n_s::Float64 = 0.0
    propellant_required_kg::Float64 = 0.0
end

Base.@kwdef struct AerobrakingControlCommand
    alpha_command::Float64 = 0.0
end

end # module CommandTypes
