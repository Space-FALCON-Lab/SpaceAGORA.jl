# import .config

function Odyssey_firing_plan(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)
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

function Earth_firing_plan(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)
    max_q = 0.6 # Max dynamic pressure
    min_q = 0.5 # Min dynamic pressure
    estimated_ρ = planet.ρ_ref * exp(-(rp-6378e3) / planet.H)
    periapsis_velocity = sqrt(planet.μ *(2/rp - 2/(ra+rp)))
    q = 0.5 * estimated_ρ * periapsis_velocity^2
    if q > max_q
        rp_new = -planet.H*log(2*min_q/(planet.ρ_ref*periapsis_velocity^2)) + 6378e3
        println("Aerobraking at ", rp, " m")
        println("New periapsis: ", rp_new, " m")
        args[:delta_v] = abs(sqrt(2*planet.μ/rp - planet.μ*2/(rp+ra)) - sqrt(2*planet.μ/rp_new - planet.μ*2/(rp_new+ra)))
        args[:phi] = deg2rad(180.0)
    else
        args[:delta_v] = 0.0
        args[:phi] = deg2rad(0.0)
    end
    return args
end

function Magellan_firing_plan(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)
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

function titan_firing_plan(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)
    titan_radius = 2575.5e3 # Radius of Titan in meters
    periapsis_radius = rp # Periapsis radius in meters
    apoapsis_radius = ra # Apoapsis radius in meters
    periapsis_altitude = periapsis_radius - titan_radius
    mu_titan = 8.981e12 # Gravitational parameter of titan in m^3/s^2
    target_periapsis_altitude = 800_000
    if periapsis_altitude < 500e3
        # Calculate delta-v required to raise periapsis to 141,000 m
        a_i = (periapsis_radius + apoapsis_radius) / 2 # Initial semi-major axis
        a_f = (target_periapsis_altitude + titan_radius + apoapsis_radius) / 2 # Final semi-major axis
        v_i = sqrt(mu_titan * (2 / apoapsis_radius - 1 / a_i))
        v_f = sqrt(mu_titan * (2 / apoapsis_radius - 1 / a_f))
        args[:delta_v] = v_f - v_i
        args[:phi] = deg2rad(180)
    else
        args[:delta_v] = 0.0
        args[:phi] = deg2rad(0)
    end
    # if args.phi == math.radians(0) and args.delta_v != 0.0:
    #     print("LOWER MANEUVER!")
    # elif args.phi == math.radians(180) and args.delta_v != 0.0:
    #     print("RAISE MANEUVER!")

    return args
end

function Venus_Express_firing_plan(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)
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