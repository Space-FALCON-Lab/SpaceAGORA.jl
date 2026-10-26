@inline function _inverse_squared_gravity_accel(pos_ii::SVector{3, Float64}, planet)::SVector{3, Float64}
    r = norm(pos_ii)
    μ = Float64(planet.μ)
    return -μ / r^2 * normalize(pos_ii)
end

@inline function _inverse_squared_j2_gravity_accel(pos_ii::SVector{3, Float64}, planet)::SVector{3, Float64}
    r = norm(pos_ii)
    μ = Float64(planet.μ)
    J2 = Float64(planet.J2)
    # J2 is normalized to the equatorial radius (matching j2_secular_rates
    # and the spherical-harmonics path), not the mean radius.
    Rp_m = Float64(planet.Rp_e)
    x, y, z = pos_ii
    r_squared = r^2
    gravity_ii_mag_spherical = -μ / r_squared
    pos_ii_hat = normalize(pos_ii)
    j2_term = SVector{3, Float64}(
        x / r * (5 * z^2 / r_squared - 1),
        y / r * (5 * z^2 / r_squared - 1),
        z / r * (5 * z^2 / r_squared - 3),
    )
    return gravity_ii_mag_spherical * pos_ii_hat + 3 / 2 * J2 * μ * Rp_m^2 / r^4 * j2_term
end
