function exponential_mars_density(config, altitude_m::Real)
    h = Float64(altitude_m)
    return config.rho_ref_kg_m3 * exp(-(h - config.h_ref_m) / config.scale_height_m)
end

function paper_heat_rate_w_cm2(config, periapsis_altitude_m::Real, periapsis_velocity_mps::Real)
    rho = max(exponential_mars_density(config, periapsis_altitude_m), eps(Float64))
    rho_nominal = max(exponential_mars_density(config, config.heat_nominal_altitude_m), eps(Float64))
    v_nominal = config.heat_nominal_velocity_mps
    heat = config.heat_nominal_w_cm2 * sqrt(rho / rho_nominal) *
           (Float64(periapsis_velocity_mps) / v_nominal)^3
    return max(0.0, heat)
end

function drag_passage_duration_s(periapsis_altitude_m::Real)
    return clamp(800.0 + 0.006 * (Float64(periapsis_altitude_m) - 92e3), 400.0, 1200.0)
end
