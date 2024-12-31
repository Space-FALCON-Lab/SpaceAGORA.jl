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

function gravity_n_bodies(et, pos_ii, p, n_body)

    primary_body_name = p.name
    n_body_name = n_body.name

    if cmp(lowercase(primary_body_name), "mars") == 0 || cmp(lowercase(primary_body_name), "jupiter") == 0 || cmp(lowercase(primary_body_name), "saturn") == 0 || cmp(lowercase(primary_body_name), "uranus") == 0 || cmp(lowercase(primary_body_name), "neptune") == 0
        primary_body_name *= "_barycenter"
    end

    if cmp(lowercase(n_body_name), "mars") == 0 || cmp(lowercase(n_body_name), "jupiter") == 0 || cmp(lowercase(n_body_name), "saturn") == 0 || cmp(lowercase(n_body_name), "uranus") == 0 || cmp(lowercase(n_body_name), "neptune") == 0
        n_body_name *= "_barycenter"
    end

    pos_primary_k = spkpos(n_body_name, et, "J2000", "none", primary_body_name)[1] * 1e3

    # println("")
    # println(pos_primary_k)
    # println("")

    pos_spacecraft_k = pos_primary_k - pos_ii
    pos_spacecraft_k_mag = norm(pos_spacecraft_k)

    g = n_body.μ * ((pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3))

    return g
end