
function def_miss(args)
    """

    """

    if args[:type_of_mission] == "Drag Passage" || args[:type_of_mission] == "Entry"
        args[:drag_passage] = 1
        args[:number_of_orbits] = 1
    elseif args[:type_of_mission] == "Orbits"
        args[:drag_passage] = 0
        args[:number_of_orbits] = args[:number_of_orbits]
    elseif args[:type_of_mission] == "Aerobraking Campaign"
        args[:drag_passage] = 0
        args[:number_of_orbits] = 1000
    end

    if args[:body_shape] == "Spacecraft"
        if args[:aerodynamic_model] == "No-Ballistic flight with axial coefficient"
            args[:aerodynamic_model] = "Mach-dependent"
            println("--AERODYNAMIC MODEL CHANGED TO: MACH-dependent - Specific for a a flat-plate--")
        end

        if args[:thermal_model] != "Maxwellian Heat Transfer"
            args[:thermal_model] = "Maxwellian Heat Transfer"
            println("--THERMAL MODEL CHANGED TO: Maxwellian Heat Transfer - Specific for a flat-plate--")
        end
    elseif args[:body_shape] == "Blunted Cone"
        if args[:aerodynamic_model] == "Mach-dependent"
            args[:aerodynamic_model] = "No-Ballistic flight with axial coefficient"
            println("--AERODYNAMIC MODEL CHANGED TO: No-Ballistic flight with axial coefficient - Specific for a blunted cone--")
        end

        if args[:thermal_model] != "Convective and Radiative"
            args[:thermal_model] = "Convective and Radiative"
            println("--THERMAL MODEL CHANGED TO: Convective and Radiative - Specific for a blunted cone--")
        end

        if args[:control_mode] != 0
            args[:control_mode] = 0
            println("--ARTICULATED SOLAR PANELS GUIDANCE NOT ALLOWED FOR BLUNTED CONE--")
        end
    end

    if args[:thrust_control] == "None"
        args[:thrust] = args[:thrust]
        args[:delta_v] = 0
    elseif args[:thrust_control] == "Aerobraking Maneuver"
        args[:thrust] = args[:thrust]
        args[:delta_v] = args[:delta_v]

        if args[:type_of_mission] == "Drag Passage"
            args[:thrust_control] = "None"
        end
    elseif args[:thrust_control] == "Drag Passage Firing"
        args[:thrust] = args[:thrust]
        args[:delta_v] = args[:delta_v]

    end

    if args[:thrust_control] != "None" && args[:thrust] == 0
        args[:thrust] = 0.1
        println("--THRUST MODIFIED TO 0.1 N--")
    end

    if Bool(args[:Odyssey_sim])
        args[:control_mode] = 0
        args[:type_of_mission] = "Aerobraking Campaign"
        args[:number_of_orbits] = 350
        args[:planet] = 1 # "Mars"
        args[:body_shape] = "Spacecraft"
        args[:dry_mass] = 411.0
        args[:prop_mass] = 50.0
        args[:α] = 90.0
        args[:inital_condition_type] = 0
        args[:thrust_control] = "Aerobraking Maneuver"

        ## For Mars Odyssey Starting at 2001-11-06
        args[:ra_initial_a] = 28559.615e3
        args[:ra_initial_b] = 30000.0e3
        args[:ra_step] = 1e12

        if args[:gravity_model] == "Inverse Squared"
            args[:hp_initial_a] = 87000 #108600
        else
            args[:hp_initial_a] = 87000 # 70000 #84200 # 86000 # 100399 # 86000 works for spherical harmonic topography (a little low, but close enough for now), 95000 for regular
        end

        args[:hp_initial_b] = 110000
        args[:hp_step] = 1e12
        args[:inclination] = 93.522
        args[:ω] = 109.7454
        args[:Ω] = 28.1517
        args[:year] = 2001
        args[:month] = 11
        args[:day] = 6

        args[:final_apoapsis] = 3390.0e3 + 503e3 # 4905.97e3  
        args[:montecarlo] = 0
        args[:drag_passage] = 0
    end

    return args
end