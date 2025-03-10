 # import .config

function Venus_Express_firing_plan(numberofpassage, args)
    if numberofpassage == 6
        args[:delta_v] = 0.428
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 36
        args[:delta_v] = 0.177
        args[:phi] = deg2rad(0)
    elseif numberofpassage == 41
        args[:delta_v] = 0.07
        args[:phi] = deg2rad(0)
    elseif numberofpassage == 45
        args[:delta_v] = 0.05
        args[:phi] = deg2rad(0)
    else
        args[:delta_v] = 0.0
        args[:phi] = deg2rad(0)
    end
    return args
end