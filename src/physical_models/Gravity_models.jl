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

    gx = (-μ * x / pos_ii_mag^3) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (1 - 5*(z/r)^2))
    gy = (-μ * y / pos_ii_mag^3) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (1 - 5*(z/r)^2))
    gz = (-μ * z / pos_ii_mag^3) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (3 - 5*(z/r)^2))

    g = [gx, gy, gz]

    return g
end