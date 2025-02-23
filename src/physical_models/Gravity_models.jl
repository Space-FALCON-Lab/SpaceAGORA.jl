using PythonCall
sys = pyimport("sys")
import .config

function gravity_const(pos_ii_mag, pos_ii, p, mass=0, vel_ii=0)
    """

    """

    μ = p.μ
    pos_ii_hat = pos_ii / pos_ii_mag

    if config.cnf.drag_state == false
        gravity_ii_mag = -μ / pos_ii_mag^2
    else
        gravity_ii_mag = -p.g_ref
    end

    g = gravity_ii_mag * pos_ii_hat

    return g
end

function gravity_invsquared(pos_ii_mag, pos_ii, p, mass=0, vel_ii=0)
    """

    """

    μ = p.μ
    pos_ii_hat = pos_ii / pos_ii_mag

    gravity_ii_mag_spherical = -μ / pos_ii_mag^2

    g = gravity_ii_mag_spherical * pos_ii_hat

    return g
end

function gravity_invsquared_J2(pos_ii_mag, pos_ii, p, mass, vel_ii=0)
    """

    """

    μ = p.μ
    J2 = p.J2

    pos_ii_hat = pos_ii ./ pos_ii_mag
    gravity_ii_mag_spherical = -μ / pos_ii_mag^2

    x = pos_ii[1]
    y = pos_ii[2]
    z = pos_ii[3]
    r = pos_ii_mag

    gx = (-μ * x / pos_ii_mag^3) * (1 + 3/2 * J2 * (p.Rp_m/r)^2 * (1 - 5*(z/r)^2))
    gy = (-μ * y / pos_ii_mag^3) * (1 + 3/2 * J2 * (p.Rp_m/r)^2 * (1 - 5*(z/r)^2))
    gz = (-μ * z / pos_ii_mag^3) * (1 + 3/2 * J2 * (p.Rp_m/r)^2 * (3 - 5*(z/r)^2))

    g = gravity_ii_mag_spherical * pos_ii_hat + 3/2 * J2 * μ * p.Rp_e^2 / r^4 * [x/r*(5*z^2/r^2 - 1), y/r*(5*z^2/r^2 - 1), z/r*(5*z^2/r^2 - 3)] 

    # g = [gx, gy, gz]

    return g
end

function gravity_GRAM(pos_ii, lat, lon, alt, p, mass, vel_ii, el_time, atmosphere, args, gram)
    """
    Function to calculate the gravity of the planet using the GRAM model
    :param pos_ii: position of the spacecraft in the inertial frame
    :param lat: latitude of the spacecraft
    :param lon: longitude of the spacecraft
    :param alt: altitude of the spacecraft
    :param timereal: current time of the simulation, used to calculate the elapsed time since the GRAM atmosphere start time
    :param atmosphere: atmosphere of the planet
    :return: gravity of the planet
    """
    # print(np.linalg.norm(pos_ii), args.EI * 1e3)
    if norm(pos_ii) - p.Rp_e > args[:EI] * 1e3
        return gravity_invsquared_J2(norm(pos_ii), pos_ii, p, mass, vel_ii)
    end
    # sys.path.append(args[:directory_Gram])
    # gram = pyimport("gram")
    position = gram.Position()
    position.lat = rad2deg(lat)
    position.lon = rad2deg(lon)
    position.height = alt*1e-3
    # position.elapsedTime = (dt.datetime(Int64(timereal.year),
    #                                         Int64(timereal.month),
    #                                         Int64(timereal.day), 
    #                                         Int64(timereal.hour),
    #                                         Int64(timereal.minute),
    #                                         timereal.second) -
    #                             dt.datetime(Int64(args.year),
    #                                         Int64(args.month),
    #                                         Int64(args.day),
    #                                         Int64(args.hours),
    #                                         Int64(args.minutes),
    #                                         args.secs)
    #                                         ).total_seconds()
    position.elapsedTime = el_time
    atmosphere.setPosition(position)
    # atmosphere.update()
    # pos = atmosphere.getPosition()
    # position.height = alt*1e-3 - pos.surfaceHeight
    # atmosphere.setPosition(position)
    atmosphere.update()
    pos = atmosphere.getPosition()

    gravity = pos.gravity
    return -pyconvert(Any, gravity) * pos_ii / norm(pos_ii)
end