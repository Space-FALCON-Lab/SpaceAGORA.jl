import .config

function titan_firing_plan(planet, apoapsis, periapsis)
    titan_radius = 2575.5e3 # Radius of Titan in meters
    periapsis_radius = periapsis # Periapsis radius in meters
    apoapsis_radius = apoapsis
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