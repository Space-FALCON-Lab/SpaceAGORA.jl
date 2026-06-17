function apoapsis_velocity(mu::Real, ra_m::Real, periapsis_altitude_m::Real, mars_radius_m::Real)
    rp = Float64(mars_radius_m) + Float64(periapsis_altitude_m)
    ra = Float64(ra_m)
    a = (ra + rp) / 2
    return sqrt(Float64(mu) * (2 / ra - 1 / a))
end

function apply_apoapsis_maneuver(config, state, action::AerobrakingAction)
    ra = state.apoapsis_radius_m
    v = apoapsis_velocity(config.mu_m3_s2, ra, state.periapsis_altitude_m, config.mars_radius_m)
    if action.lowers_periapsis
        v -= action.magnitude_mps
    elseif action.raises_periapsis
        v += action.magnitude_mps
    end
    energy = 0.5 * v^2 - config.mu_m3_s2 / ra
    semi_major = -config.mu_m3_s2 / (2 * energy)
    periapsis_altitude = 2 * semi_major - ra - config.mars_radius_m
    return periapsis_altitude
end
