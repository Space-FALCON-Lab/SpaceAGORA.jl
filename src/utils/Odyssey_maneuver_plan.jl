import .config


function Odyssey_firing_plan(numberofpassage, args)
    if numberofpassage == 7
        args[:delta_v] = 0.14
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 14
        args[:delta_v] = 0.15
        args[:phi] = 0
    elseif numberofpassage == 26
        args[:delta_v] = 0.1
        args[:phi] = 0
    elseif numberofpassage == 30
        args[:delta_v] = 0.1
        args[:phi] = 0
    elseif numberofpassage == 35
        args[:delta_v] = 0.2
        args[:phi] = 0
    elseif numberofpassage == 47
        args[:delta_v] = 0.2
        args[:phi] = 0
    elseif numberofpassage == 55
        args[:delta_v] = 0.3
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 69
        args[:delta_v] = 0.15
        args[:phi] = 0
    elseif numberofpassage == 72
        args[:delta_v] = 0.15
        args[:phi] = 0
    elseif numberofpassage == 80
        args[:delta_v] = 0.15
        args[:phi] = 0
    elseif numberofpassage == 87
        args[:delta_v] = 0.14 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 110
        args[:delta_v] = 0.14 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 127
        args[:delta_v] = 1.0 # /3
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 160
        args[:delta_v] = 0.84 #/2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 178
        args[:delta_v] = 0.6  #/2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 194
        args[:delta_v] = 0.84 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 210
        args[:delta_v] = 0.6 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 222
        args[:delta_v] = 0.6 #/2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 238
        args[:delta_v] = 1.2 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 251
        args[:delta_v] = 1.0 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 263
        args[:delta_v] = 1.0 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 274
        args[:delta_v] = 1.2 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 287
        args[:delta_v] = 1.0 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 299
        args[:delta_v] = 1.0 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 311
        args[:delta_v] = 1.2 # /2
        args[:phi] = deg2rad(180)
    else
        args[:delta_v] = 0.0
        args[:phi] = 0.0
    end

    return args
end