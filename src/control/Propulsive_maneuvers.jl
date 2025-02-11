include("../utils/Odyssey_maneuver_plan.jl")

import .config

function no_maneuver(t0, thrust_mag, Δv, args, index_phase_aerobraking)
    thrust = 0
    return thrust
end

function abms(t0, thrust_mag, Δv, args, index_phase_aerobraking)
    if index_phase_aerobraking == 0
        thrust = thrust_mag
    else
        thrust = 0
    end

    return thrust
end

function deceleration_drag_passage(t0, thrust_mag, Δv, args, index_phase_aerobraking)
    if config.cnf.drag_state == true
        thrust = thrust_mag
    else
        thrust = 0
    end

    return thrust
end