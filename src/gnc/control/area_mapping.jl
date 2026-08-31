#=
"""
    Exposed area to AOA relation.

    The controller variable is exposed area. The panel command is recovered
    from A(theta) = S_bus + S_SP sin(theta), with theta in [0, pi/2].
"""
=#
function commanded_area_fraction(config::AerobrakingMPCConfig, area_m2::Real)
    denom = max(config.controllable_area_m2, eps(Float64))
    return clamp((Float64(area_m2) - config.bus_reference_area_m2) / denom, 0.0, 1.0)
end

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

function commanded_area_from_alpha(
    config::AerobrakingMPCConfig,
    alpha_rad::Real;
    min_alpha_rad::Real,
    max_alpha_rad::Real,
)
    alpha = clamp(Float64(alpha_rad), Float64(min_alpha_rad), Float64(max_alpha_rad))
    return config.bus_reference_area_m2 + config.controllable_area_m2 * sin(alpha)
end

function apply_commanded_area!(
    spacecraft,
    config::AerobrakingMPCConfig,
    area_m2::Real;
    controlled_panel_links,
    min_alpha_rad::Real,
    max_alpha_rad::Real,
)
    α = alpha_from_commanded_area(
        config,
        area_m2;
        min_alpha_rad=min_alpha_rad,
        max_alpha_rad=max_alpha_rad,
    )
    links = _maybe_property(spacecraft, :links, nothing)
    links === nothing && throw(ArgumentError("Cannot set panel angles on spacecraft without links."))
    for link_index in controlled_panel_links
        links[Int(link_index)].α = α
    end
    return α
end
