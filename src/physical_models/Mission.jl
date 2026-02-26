@inline _legacy_key(x) = lowercase(strip(replace(string(x), r"[_-]+" => " ")))

if !isdefined(@__MODULE__, :_deprecated_mission_dict_input_warned)
    const _deprecated_mission_dict_input_warned = Ref(Set{Symbol}())
end

@inline _warn_deprecated_legacy_dict_enabled() = get(ENV, "SPACEAGORA_WARN_DEPRECATED_CONFIG", "1") == "1"

function _warn_deprecated_mission_dict_value!(key::Symbol, value)
    if !_warn_deprecated_legacy_dict_enabled() || (key in _deprecated_mission_dict_input_warned[])
        return nothing
    end
    push!(_deprecated_mission_dict_input_warned[], key)
    @warn "mission_def: legacy string/symbol value for :$key ($(repr(value))) is deprecated; migrate to typed config API."
    return nothing
end

@inline function _legacy_norm_mission_value(mission::Dict{Symbol, Any}, key::Symbol, default)
    has_user_value = haskey(mission, key)
    value = has_user_value ? mission[key] : default
    if has_user_value && (value isa AbstractString || value isa Symbol)
        _warn_deprecated_mission_dict_value!(key, value)
    end
    return _legacy_key(value)
end

@inline function _legacy_mission_planet_id(x)::Int
    key = _legacy_key(x)
    if x isa Integer
        pid = Int(x)
        if pid in (0, 1, 2, 3, 7)
            return pid
        end
    elseif key == "earth"
        return 0
    elseif key == "mars"
        return 1
    elseif key == "venus"
        return 2
    elseif key == "sun"
        return 3
    elseif key == "titan"
        return 7
    end
    return 1 # legacy default (Mars)
end

function mission_def(mission::Dict{Symbol, Any})

    e, d, l, a = 0, 0, 0, 1     # e = Entry, d = Descent, l = Landing, a = Aerobraking : 0 - No, 1 - Yes
    planet_val = get(mission, :Planet, 1)
    if haskey(mission, :Planet) && (planet_val isa AbstractString || planet_val isa Symbol)
        _warn_deprecated_mission_dict_value!(:Planet, planet_val)
    end
    p = _legacy_mission_planet_id(planet_val)

    grav_key = _legacy_norm_mission_value(mission, :Gravity_Model, "Inverse Squared")
    dens_key = _legacy_norm_mission_value(mission, :Density_Model, "Exponential")
    aero_key = _legacy_norm_mission_value(mission, :Aerodynamic_Model, "Cd and Cl Constant")
    shape_key = _legacy_norm_mission_value(mission, :Shape, "Spacecraft")
    thermal_key = _legacy_norm_mission_value(mission, :Thermal_Model, "Convective and Radiative")
    firing_key = _legacy_norm_mission_value(mission, :Firings, "None")

    M = Mission(e, d, l, a, p)

    # Gravity Model Selection
    if grav_key == "constant"
        gm = LegacyGravityConstant
    elseif grav_key == "inverse squared"
        gm = LegacyGravityInverseSquared
    elseif grav_key == "inverse squared and j2 effect"
        gm = LegacyGravityInverseSquaredJ2
    elseif grav_key == "gram"
        gm = LegacyGravityGRAM
    else
        gm = LegacyGravityInverseSquared
    end

    # Density Model Selection
    if dens_key == "constant"
        dm = LegacyDensityConstant
    elseif dens_key == "exponential"
        dm = LegacyDensityExponential
    elseif dens_key == "no density"
        dm = LegacyDensityNoDensity
    elseif dens_key == "gram"
        dm = LegacyDensityGRAM
    elseif dens_key == "nrlmsise"
        dm = LegacyDensityNRLMSISE
    else
        dm = LegacyDensityExponential
    end

    # Wind
    wm = Int64(get(mission, :Wind, 0))

    # Aerodynamic Model Selection
    if aero_key in ("cd and cl constant",)
        am = LegacyAeroCdClConstant
    elseif aero_key == "diffusive" || (aero_key == "mach dependent" && shape_key == "spacecraft")
        am = LegacyAeroDiffusive
    elseif aero_key in ("no balistic flight with axial coefficent", "no ballistic flight with axial coefficient") &&
           shape_key == "blunted cone"
        am = LegacyAeroNoBallisticAxial
    else
        am = LegacyAeroCdClConstant
    end
    
    # Control Mode
    control_val = get(mission, :Control, 0)
    if control_val == 3
        cm = 3
    elseif control_val == 2
        cm = 2
    elseif control_val == 1
        cm = 1
    else
        cm = 0
    end

    # Thrust Control
    if firing_key == "none"
        tc = LegacyThrustNone
    elseif firing_key == "aerobraking maneuver"
        tc = LegacyThrustAerobrakingManeuver
    elseif firing_key == "drag passage firing"
        tc = LegacyThrustDragPassageFiring
    else
        tc = LegacyThrustNone
    end

    # Thermal Model
    if thermal_key == "convective and radiative"
        tm = LegacyThermalConvectiveRadiative
    elseif thermal_key in ("maxwellian heat transfer", "shaaf and chambre")
        tm = LegacyThermalMaxwellian
    else
        tm = LegacyThermalConvectiveRadiative
    end

    # MonteCarlo
    mc = Int64(get(mission, :Monte_Carlo, 0))


    ip = InitialParameters(M, gm, dm, wm, am, tm, cm, tc, mc)

    return ip
end
