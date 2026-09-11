module LaserTerminal

export LaserTerminalActuator, LaserTerminalCommand, apply_laser_terminal_command!

Base.@kwdef mutable struct LaserTerminalActuator
    enabled::Bool = true
    boresight_body::NTuple{3, Float64} = (1.0, 0.0, 0.0)
    max_pointing_error_rad::Float64 = 0.0
end

Base.@kwdef struct LaserTerminalCommand
    target_id::String = ""
    tx_power_w::Float64 = 0.0
    acquire_link::Bool = false
end

function apply_laser_terminal_command!(::LaserTerminalActuator, ::LaserTerminalCommand, ::Float64)
    throw(ErrorException("Not implemented: apply_laser_terminal_command!(::LaserTerminalActuator, ::LaserTerminalCommand, ::Float64)"))
end

end # module LaserTerminal
