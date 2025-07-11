
using PythonCall

sys = pyimport("sys")

function gravity_const(pos_ii_mag::Float64, pos_ii::SVector{3, Float64}, p)
    """

    """

    μ = p.μ
    pos_ii_hat = normalize(pos_ii)

    if config.cnf.drag_state == false
        gravity_ii_mag = -μ / pos_ii_mag^2
    else
        gravity_ii_mag = -p.g_ref
    end

    return gravity_ii_mag * pos_ii_hat
end

function gravity_invsquared(pos_ii_mag::Float64, pos_ii::SVector{3, Float64}, p)
    """

    """
    return -p.μ / pos_ii_mag^2 * normalize(pos_ii)
end

function gravity_invsquared_J2(r::Float64, pos_ii::SVector{3, Float64}, p)
    """

    """

    μ = p.μ
    J2 = p.J2

    pos_ii_hat = normalize(pos_ii)
    r_squared = r^2
    gravity_ii_mag_spherical = -μ / r_squared

    x,y,z = pos_ii

    # gx = (-μ * x / r_cubed) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (1 - 5*(z/r)^2))
    # gy = (-μ * y / r_cubed) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (1 - 5*(z/r)^2))
    # gz = (-μ * z / r_cubed) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (3 - 5*(z/r)^2))

    return gravity_ii_mag_spherical * pos_ii_hat + 3/2 * J2 * μ * p.Rp_m^2 / r^4 * [x/r*(5*z^2/r_squared - 1), y/r*(5*z^2/r_squared - 1), z/r*(5*z^2/r_squared - 3)] 
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
    if norm(pos_ii) - p.Rp_e > args[:EI] * 1e3
        return gravity_invsquared_J2(norm(pos_ii), pos_ii, p)
    end

    position = gram.Position()
    position.lat = rad2deg(lat)
    position.lon = rad2deg(lon)
    position.height = alt*1e-3

    position.elapsedTime = el_time
    atmosphere.setPosition(position)
    
    atmosphere.update()
    pos = atmosphere.getPosition()

    gravity = pos.gravity
    return -pyconvert(Any, gravity) * pos_ii / norm(pos_ii)
end