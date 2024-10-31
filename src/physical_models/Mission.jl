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
    p = 1                       # 0 - Earth, 1 - Mars, 2 - Venus

    M = Mission(e, d, l, a, p)

    # Gravity Model Selection
    if mission[:Gravity_Model] == "constant" || mission[:Gravity_Model] == "Constant"
        gm = 0
    elseif mission[:Gravity_Model] == "Inverse squared" || mission[:Gravity_Model] == "inverse squared" || mission[:Gravity_Model] == "Inverse Squared"
        gm = 1
    elseif mission[:Gravity_Model] == "Inverse Squared and J2 effect" || mission[:Gravity_Model] == "inverse squared and J2 effect" || mission[:Gravity_Model] == "Inverse quared and J2 effect"
        gm = 2
    else
        gm = 1
    end

    # Density Model Selection
    if mission[:Density_Model] == "constant" || mission[:Density_Model] == "Constant"
        dm = 0
    elseif mission[:Density_Model] == "expopnential" || mission[:Density_Model] == "Exponential"
        dm = 1
    elseif mission[:Density_Model] == "No-Density" || mission[:Density_Model] == "No-density"
        dm = 2
    elseif mission[:Density_Model] == "MARSGram" || mission[:Density_Model] == "MarsGram"
        dm = 3
    else
        dm = 1
    end

    # Wind
    wm = Int64(mission[:Wind])

    # Aerodynamic Model Selection
    if mission[:Aerodynamic_Model] == "Cd and Cl Constant" || mission[:Aerodynamic_Model] == "Cd and Cl constant"
        am = 0
    elseif mission[:Aerodynamic_Model] == "Diffusive" || mission[:Aerodynamic_Model] == "Mach-dependent" && mission[:Aerodynamic_Model] == "Diffusive" || mission[:Aerodynamic_Model] == "Spacecraft"
        am = 1
    elseif mission[:Aerodynamic_Model] == "No-Balistic flight with axial coefficent" || mission[:Aerodynamic_Model] == "No-ballistic flight with axial coefficient" && mission[:Shape] == "Blunted Cone"
        am = 2
    else
        am = 0
    end
    
    # Control Mode
    if mission[:Control] == 3
        cm = 1
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
    if mission[:Thermal_Model] == "convective and radiative" || mission[:Thermal_model] == "Convective and Radiative"
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