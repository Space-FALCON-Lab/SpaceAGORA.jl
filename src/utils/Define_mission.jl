
@inline _legacy_token(x) = lowercase(strip(replace(string(x), r"[_-]+" => " ")))

@inline function _legacy_mission_kind(x)::Symbol
    key = _legacy_token(x)
    if key in ("drag passage",)
        return :drag_passage
    elseif key in ("entry",)
        return :entry
    elseif key in ("orbits", "orbit", "missionorbits")
        return :orbits
    elseif key in ("time", "missiontime")
        return :time
    elseif key in ("aerobraking campaign",)
        return :campaign
    end
    return :unknown
end

function def_miss(args)
    """

    """

    mission_kind = _legacy_mission_kind(get(args, :type_of_mission, "Time"))
    if mission_kind == :drag_passage || mission_kind == :entry
        args[:drag_passage] = 1
        args[:number_of_orbits] = 1
    elseif mission_kind == :orbits || mission_kind == :time
        args[:drag_passage] = 0
        # args[:number_of_orbits] = args[:number_of_orbits]
    elseif mission_kind == :campaign
        args[:drag_passage] = 0
        args[:number_of_orbits] = 1000
    end

    body_shape = _legacy_token(get(args, :body_shape, "Spacecraft"))
    aerodynamic_model = _legacy_token(get(args, :aerodynamic_model, "Mach-dependent"))
    thermal_model = _legacy_token(get(args, :thermal_model, "Maxwellian Heat Transfer"))
    thrust_control = _legacy_token(get(args, :thrust_control, "None"))

    if body_shape == "spacecraft"
        if aerodynamic_model in ("no ballistic flight with axial coefficient", "no balistic flight with axial coefficent")
            args[:aerodynamic_model] = "Mach-dependent"
            println("--AERODYNAMIC MODEL CHANGED TO: MACH-dependent - Specific for a flat-plate--")
        end

        if thermal_model != "maxwellian heat transfer"
            args[:thermal_model] = "Maxwellian Heat Transfer"
            println("--THERMAL MODEL CHANGED TO: Maxwellian Heat Transfer - Specific for a flat-plate--")
        end
    elseif body_shape == "blunted cone"
        if aerodynamic_model == "mach dependent"
            args[:aerodynamic_model] = "No-Ballistic flight with axial coefficient"
            println("--AERODYNAMIC MODEL CHANGED TO: No-Ballistic flight with axial coefficient - Specific for a blunted cone--")
        end

        if thermal_model != "convective and radiative"
            args[:thermal_model] = "Convective and Radiative"
            println("--THERMAL MODEL CHANGED TO: Convective and Radiative - Specific for a blunted cone--")
        end

        if get(args, :control_mode, 0) != 0
            args[:control_mode] = 0
            println("--ARTICULATED SOLAR PANELS GUIDANCE NOT ALLOWED FOR BLUNTED CONE--")
        end
    end

    if thrust_control == "none"
        args[:thrust] = args[:thrust]
        args[:delta_v] = 0
    elseif thrust_control == "aerobraking maneuver"
        args[:thrust] = args[:thrust]
        args[:delta_v] = args[:delta_v]

        if mission_kind == :drag_passage
            args[:thrust_control] = "None"
        end
    elseif thrust_control == "drag passage firing"
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
        # args[:dry_mass] = 411.0
        # args[:prop_mass] = 50.0
        args[:α] = 90.0
        args[:inital_condition_type] = 0
        args[:thrust_control] = "Aerobraking Maneuver"

        ## For Mars Odyssey Starting at 2001-11-06
        args[:ra_initial_a] = 28559.615e3 # 28814.747e3
        args[:ra_initial_b] = 30000.0e3
        args[:ra_step] = 1e12

        if args[:gravity_model] == "Inverse Squared"
            args[:hp_initial_a] = 87000 #108600
        else
            args[:hp_initial_a] = 87000 #88200 #88500 # 70000 #84200 # 86000 # 100399 # 86000 works for spherical harmonic topography (a little low, but close enough for now), 95000 for regular
        end

        args[:hp_initial_b] = 110000
        args[:hp_step] = 1e12
        args[:inclination] = 93.522
        args[:ω] = 109.7454
        args[:Ω] = 28.1517
        # args[:inclination] = 93.50725189698
        # args[:ω] = 109.926987321
        # args[:Ω] = 28.123464191
        args[:year] = 2001
        args[:month] = 11
        args[:day] = 6

        args[:final_apoapsis] = 3390.0e3 + 503e3 # 4905.97e3  
        args[:montecarlo] = 0
        args[:drag_passage] = 0
    end

    return args
end
