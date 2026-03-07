module CommandTypes

export PropulsiveManeuverCommand, AerobrakingControlCommand

Base.@kwdef struct PropulsiveManeuverCommand
    valid::Bool = false
    delta_v_mps::Float64 = 0.0
    direction_rad::Float64 = 0.0
    source_orbit::Int64 = 0
end

Base.@kwdef struct AerobrakingControlCommand
    alpha_command::Float64 = 0.0
end

end # module CommandTypes
