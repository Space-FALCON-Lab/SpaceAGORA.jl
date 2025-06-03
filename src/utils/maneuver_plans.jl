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
    elseif numberofpassage == 54
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
        args[:delta_v] = 0.15 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 110
        args[:delta_v] = 0.15 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 128
        args[:delta_v] = 1.0 # /3
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 161
        args[:delta_v] = 0.84 #/2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 179
        args[:delta_v] = 0.6  #/2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 195
        args[:delta_v] = 0.84 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 211
        args[:delta_v] = 0.6 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 223
        args[:delta_v] = 0.6 #/2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 239
        args[:delta_v] = 1.2 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 251
        args[:delta_v] = 1.0 # /2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 263
        args[:delta_v] = 1.2 # /2
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

# Use when starting from true beginning
function Odyssey_firing_plan_true_beginning(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)
    if numberofpassage == 7
        args[:delta_v] = 0.549
        args[:phi] = 0.0
    elseif numberofpassage == 10
        args[:delta_v] = 0.549
        args[:phi] = 0.0
    elseif numberofpassage == 12
        args[:delta_v] = 0.301
        args[:phi] = 0.0
    elseif numberofpassage == 14
        args[:delta_v] = 0.151
        args[:phi] = 0.0
    elseif numberofpassage == 16
        args[:delta_v] = 0.303
        args[:phi] = 0.0
    elseif numberofpassage == 25
        args[:delta_v] = 0.145
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 32
        args[:delta_v] = 0.151
        args[:phi] = 0.0
    elseif numberofpassage == 45
        args[:delta_v] = 0.099
        args[:phi] = 0.0
    elseif numberofpassage == 49
        args[:delta_v] = 0.099
        args[:phi] = 0.0
    elseif numberofpassage == 53
        args[:delta_v] = 0.2
        args[:phi] = 0.0
    elseif numberofpassage == 66
        args[:delta_v] = 0.2
        args[:phi] = 0.0
    elseif numberofpassage == 74
        args[:delta_v] = 0.297
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 88
        args[:delta_v] = 0.149
        args[:phi] = 0.0
    elseif numberofpassage == 91
        args[:delta_v] = 0.149
        args[:phi] = 0.0
    elseif numberofpassage == 99
        args[:delta_v] = 0.147
        args[:phi] = 0.0
    elseif numberofpassage == 105
        args[:delta_v] = 0.144
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 129
        args[:delta_v] = 0.144
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 146
        args[:delta_v] = 1.0
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 179
        args[:delta_v] = 0.844
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 197
        args[:delta_v] = 0.6
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 213
        args[:delta_v] = 0.844
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 229
        args[:delta_v] = 0.6
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 241
        args[:delta_v] = 0.6
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 257
        args[:delta_v] = 1.2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 270
        args[:delta_v] = 1.0
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 282
        args[:delta_v] = 1.2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 295
        args[:delta_v] = 1.2
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 306
        args[:delta_v] = 1.0
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 318
        args[:delta_v] = 1.0
        args[:phi] = deg2rad(180)
    elseif numberofpassage == 330
        args[:delta_v] = 1.2
        args[:phi] = deg2rad(180)
    else
        args[:delta_v] = 0.0
        args[:phi] = 0.0
    end

    return args
end
function Earth_firing_plan(planet=nothing, ra=0.0, rp=0.0, numberofpassage=0.0, args=nothing)
    # max_q = 100000 # Max heat rate
    # min_q = 99000 # Min heat rate
    # estimated_ρ = planet.ρ_ref * exp(-(rp-6378e3) / planet.H)
    # periapsis_velocity = sqrt(planet.μ *(2/rp - 2/(ra+rp)))
    # q = 0.5 * estimated_ρ * periapsis_velocity^3
    # println("Estimated q: ", q, " W/cm^2")
    if rp < 120e3 + planet.Rp_e
        # rp_new = -planet.H*log(2*min_q/(planet.ρ_ref*periapsis_velocity^3)) + 6378e3
        # println("Aerobraking at ", rp, " m")
        # println("New periapsis: ", rp_new, " m")
        target_periapsis_altitude = 140e3 # Target periapsis altitude in meters
        a_i = (rp + ra) / 2 # Initial semi-major axis
        a_f = (target_periapsis_altitude + planet.Rp_e + ra) / 2 # Final semi-major axis
        v_i = sqrt(planet.μ * (2 / ra - 1 / a_i))
        v_f = sqrt(planet.μ * (2 / ra - 1 / a_f))
        args[:delta_v] = v_f - v_i
        # args[:delta_v] = abs(sqrt(2*planet.μ/rp - planet.μ*2/(rp+ra)) - sqrt(2*planet.μ/rp_new - planet.μ*2/(rp_new+ra)))
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
        args[:phi] = 0.0
    elseif numberofpassage == 508
        args[:delta_v] = one_n
        args[:phi] = 0.0
    elseif numberofpassage == 599
        args[:delta_v] = one_n
        args[:phi] = 0.0
    else
        args[:delta_v] = 0.0
        args[:phi] = 0.0
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
    target_min_periapsis_altitude = 550_000
    if periapsis_altitude < target_min_periapsis_altitude || abs(periapsis_altitude - target_min_periapsis_altitude) < 5e3
        # Calculate delta-v required to raise periapsis to 141,000 m
        a_i = (periapsis_radius + apoapsis_radius) / 2 # Initial semi-major axis
        a_f = (target_periapsis_altitude + titan_radius + apoapsis_radius) / 2 # Final semi-major axis
        v_i = sqrt(mu_titan * (2 / apoapsis_radius - 1 / a_i))
        v_f = sqrt(mu_titan * (2 / apoapsis_radius - 1 / a_f))
        args[:delta_v] = v_f - v_i
        args[:phi] = deg2rad(180)
    else
        args[:delta_v] = 0.0
        args[:phi] = 0.0
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
        args[:phi] = 0.0
    elseif numberofpassage == 41
        args[:delta_v] = 0.07
        args[:phi] = 0.0
    elseif numberofpassage == 45
        args[:delta_v] = 0.05
        args[:phi] = 0.0
    else
        args[:delta_v] = 0.0
        args[:phi] = 0.0
    end
    return args
end