import .config

function Magellan_firing_plan(numberofpassage, args)
    one_n = 0.34 # Δv of a 1-n maneuver
    two_n = 0.68 # Δv of a 2-n maneuver
    half_n = 0.17 # Δv of a 1/2-n maneuver

    if numberofpassage == 50
        args[:delta_v] = half_n
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 147
        args[:delta_v] = half_n
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 185
        args[:delta_v] = half_n
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 212
        args[:delta_v] = half_n
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 238
        args[:delta_v] = half_n
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 297
        args[:delta_v] = one_n
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 444
        args[:delta_v] = half_n
        args[:phi] = deg2rad(0)
    elseif numberofpassage == 508
        args[:delta_v] = one_n
        args[:phi] = deg2rad(0)
    elseif numberofpassage == 599
        args[:delta_v] = one_n
        args[:phi] = deg2rad(0)
    else
        args[:delta_v] = 0.0
        args[:phi] = deg2rad(0)
    end

    return args
end