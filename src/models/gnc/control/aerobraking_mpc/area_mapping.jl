# Map exposed area to panel angle using A(θ) = S_bus + S_SP sin(θ).
function commanded_area_fraction(config::AerobrakingMPCConfig, area_m2::Real)
    denom = max(config.controllable_area_m2, eps(Float64))
    return clamp((Float64(area_m2) - config.bus_reference_area_m2) / denom, 0.0, 1.0)
end

"""
    alpha_from_commanded_area(config, area_m2; min_alpha_rad, max_alpha_rad) -> Float64

Panel angle of attack (rad) giving `area_m2 = bus + controllable * sin(α)`, clamped to the range.
"""
function alpha_from_commanded_area(
    config::AerobrakingMPCConfig,
    area_m2::Real;
    min_alpha_rad::Real,
    max_alpha_rad::Real,
)
    exposed_fraction = commanded_area_fraction(config, area_m2)
    alpha = asin(exposed_fraction)
    return clamp(alpha, Float64(min_alpha_rad), Float64(max_alpha_rad))
end

"""
    commanded_area_from_alpha(config, alpha_rad; min_alpha_rad, max_alpha_rad) -> Float64

Exposed area (m²) `bus + controllable * sin(α)` for a panel angle of attack.
"""
function commanded_area_from_alpha(
    config::AerobrakingMPCConfig,
    alpha_rad::Real;
    min_alpha_rad::Real,
    max_alpha_rad::Real,
)
    alpha = clamp(Float64(alpha_rad), Float64(min_alpha_rad), Float64(max_alpha_rad))
    return config.bus_reference_area_m2 + config.controllable_area_m2 * sin(alpha)
end
