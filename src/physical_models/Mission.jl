mutable struct Mission
    e::Int64
    d::Int64
    l::Int64
    a::Int64
    planet::Int64
end

mutable struct InitalParameters
    M::Mission
    gm::Int64
    dm::Int64
    wm::Int64
    am::Int64
    tm::Int64
    cm::Int64
    tc::Int64
    mc::Int64
end

function mission_def(mission)

    e, d, l, a = 0, 0, 0, 1     # e = Entry, d = Descent, l = Landing, a = Aerobraking : 0 - No, 1 - Yes

    if (mission[:Planet] == 0 || (typeof(mission[:Planet]) == String && cmp(lowercase(mission[:Planet]), "earth") == 0)) # Earth
        p = 0
    elseif (mission[:Planet] == 1 || (typeof(mission[:Planet]) == String && cmp(lowercase(mission[:Planet]), "mars") == 0)) # Mars
        p = 1
    elseif (mission[:Planet] == 2 || (typeof(mission[:Planet]) == String && cmp(lowercase(mission[:Planet]), "venus") == 0)) # Venus
        p = 2
    elseif (mission[:Planet] == 3 || (typeof(mission[:Planet]) == String && cmp(lowercase(mission[:Planet]), "sun") == 0)) # Sun
        p = 3
    elseif ((mission[:Planet] == 7 || (typeof(mission[:Planet]) == String && cmp(lowercase(mission[:Planet]), "titan") == 0))) # Titan
        p = 7
    else
        p = 1
    end

    M = Mission(e, d, l, a, p)

    # Gravity Model Selection
    if uppercase(mission[:Gravity_Model]) == "CONSTANT"
        gm = 0
    elseif uppercase(mission[:Gravity_Model]) == "INVERSE SQUARED"
        gm = 1
    elseif uppercase(mission[:Gravity_Model]) == "INVERSE SQUARED AND J2 EFFECT"
        gm = 2
    elseif uppercase(mission[:Gravity_Model]) == "GRAM"
        gm = 3
    else
        gm = 1
    end

    # Density Model Selection
    if uppercase(mission[:Density_Model]) == "CONSTANT"
        dm = 0
    elseif uppercase(mission[:Density_Model]) == "EXPONENTIAL"
        dm = 1
    elseif uppercase(mission[:Density_Model]) == "NO-DENSITY"
        dm = 2
    elseif uppercase(mission[:Density_Model]) == "GRAM"
        dm = 3
    elseif uppercase(mission[:Density_Model]) == "NRLMSISE"
        dm = 4
    else
        dm = 1
    end

    # Wind
    wm = Int64(mission[:Wind])

    # Aerodynamic Model Selection
    if mission[:Aerodynamic_Model] == "Cd and Cl Constant" || mission[:Aerodynamic_Model] == "Cd and Cl constant"
        am = 0
    elseif mission[:Aerodynamic_Model] == "Diffusive" || mission[:Aerodynamic_Model] == "Mach-dependent" &&  mission[:Shape] == "Spacecraft"
        am = 1
    elseif mission[:Aerodynamic_Model] == "No-Balistic flight with axial coefficent" || mission[:Aerodynamic_Model] == "No-ballistic flight with axial coefficient" && mission[:Shape] == "Blunted Cone"
        am = 2
    else
        am = 0
    end
    
    # Control Mode
    if mission[:Control] == 3
        cm = 3
    elseif mission[:Control] == 2
        cm = 2
    elseif mission[:Control] == 1
        cm = 1
    else
        cm = 0
    end

    # Thrust Control
    if mission[:Firings] == "None"
        tc = 0
    elseif mission[:Firings] == "Aerobraking Maneuver"
        tc = 1
    elseif mission[:Firings] == "Drag Passage Firing"
        tc = 2
    else
        tc = 0
    end

    # Thermal Model
    if mission[:Thermal_Model] == "convective and radiative" || mission[:Thermal_Model] == "Convective and Radiative"
        tm = 1
    elseif mission[:Thermal_Model] == "Maxwellian Heat Transfer" || mission[:Thermal_Model] == "Shaaf and Chambre"
        tm = 2
    else
        tm = 1
    end

    # MonteCarlo
    mc = Int64(mission[:Monte_Carlo])


    ip = InitalParameters(M, gm, dm, wm, am, tm, cm, tc, mc)

    return ip
end