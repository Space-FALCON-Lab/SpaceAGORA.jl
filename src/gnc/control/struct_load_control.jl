using Roots

@inline function _edg_controlled_panel_area(
    spacecraft,
    controlled_panel_links::Tuple{Vararg{Int}},
)::Float64
    area = 0.0
    for idx in controlled_panel_links
        if 1 <= idx <= length(spacecraft.links)
            area += max(0.0, Float64(spacecraft.links[idx].ref_area))
        end
    end
    return area
end

@inline function _edg_legacy_plate_coefficients(speed_ratio::Float64, alpha::Float64, sigma::Float64)
    s = max(speed_ratio, eps(Float64))
    s_sin = s * sin(alpha)
    common = exp(-(s_sin)^2)
    erf_common = 1.0 + erf(s_sin)
    cn = (
        (((2.0 - sigma) / sqrt(pi)) * s_sin + sigma / 2.0) * common +
        ((2.0 - sigma) * (s_sin^2 + 0.5) + sigma * sqrt(pi) * s_sin / 2.0) * erf_common
    ) / s^2
    ca = sigma * cos(alpha) * (common + sqrt(pi) * s_sin * erf_common) / (sqrt(pi) * s)
    cl = cn * cos(alpha) - ca * sin(alpha)
    cd = ca * cos(alpha) + cn * sin(alpha)
    return cl, cd
end


function _edg_legacy_spacecraft_aero_coefficients(
    spacecraft,
    speed_ratio::Float64,
    controlled_panel_links::Tuple{Vararg{Int}},
    controlled_alpha::Float64,
)
    panel_area = _edg_controlled_panel_area(spacecraft, controlled_panel_links)
    total_area = sum(max(0.0, Float64(link.ref_area)) for link in spacecraft.links)
    bus_area = max(0.0, total_area - panel_area)
    sigma = Float64(spacecraft.root.reflection_coefficient)
    cl_panel, cd_panel = _edg_legacy_plate_coefficients(speed_ratio, controlled_alpha, sigma)
    cl_bus, cd_bus = _edg_legacy_plate_coefficients(speed_ratio, pi / 2, sigma)
    inv_area = 1.0 / max(total_area, eps(Float64))
    return (
        (cl_panel * panel_area + cl_bus * bus_area) * inv_area,
        (cd_panel * panel_area + cd_bus * bus_area) * inv_area,
    )
end

function _edg_dynamic_spacecraft_aero_coefficients(
    spacecraft,
    temperature::Float64,
    speed_ratio::Float64,
    controlled_panel_links::Tuple{Vararg{Int}},
    controlled_alpha::Float64,
)
    total_area = 0.0
    cl_area = 0.0
    cd_area = 0.0
    for (idx, link) in pairs(spacecraft.links)
        area = max(0.0, Float64(link.ref_area))
        alpha = idx in controlled_panel_links ? controlled_alpha : pi / 2
        cl, cd, _, _, _, _ = aerodynamic_coefficient_fM(
            link,
            temperature,
            speed_ratio,
            alpha,
            0.0,
            0.0,
        )
        total_area += area
        cl_area += cl * area
        cd_area += cd * area
    end
    inv_area = 1.0 / max(total_area, eps(Float64))
    return cl_area * inv_area, cd_area * inv_area
end


function _energy_depletion_struct_drag_area(
    spacecraft,
    _temperature::Float64,
    speed_ratio::Float64,
    controlled_panel_links::Tuple{Vararg{Int}},
    controlled_alpha::Float64,
    _config,
)::Float64
    total_area = sum(max(0.0, Float64(link.ref_area)) for link in spacecraft.links)
    _, cd = _edg_legacy_spacecraft_aero_coefficients(
        spacecraft,
        speed_ratio,
        controlled_panel_links,
        controlled_alpha,
    )
    return cd * total_area
end

function _energy_depletion_struct_load_root_alpha(
    config,
    env,
    spacecraft,
    controlled_panel_links::Tuple{Vararg{Int}},
    base_alpha::Float64,
)::Float64
    limit = config.structural_load_limit_pa
    q = env.dynamic_pressure
    if !(isfinite(limit) && limit > 0.0 && isfinite(q) && q > 0.0)
        return base_alpha
    end

    min_alpha = _edg_constraint_min_alpha(config)
    max_alpha = clamp(base_alpha, min_alpha, config.max_alpha_rad)
    temperature = env.temperature
    speed_ratio = max(env.molecular_speed_ratio, eps(Float64))

    reference_drag_area = _energy_depletion_struct_drag_area(
        spacecraft,
        temperature,
        speed_ratio,
        controlled_panel_links,
        config.max_alpha_rad,
        config,
    )
    drag_limit = limit * reference_drag_area

    drag_at_max = q * _energy_depletion_struct_drag_area(
        spacecraft,
        temperature,
        speed_ratio,
        controlled_panel_links,
        max_alpha,
        config,
    )
    drag_at_min = q * _energy_depletion_struct_drag_area(
        spacecraft,
        temperature,
        speed_ratio,
        controlled_panel_links,
        min_alpha,
        config,
    )

    if drag_at_max <= drag_limit
        return max_alpha
    elseif drag_at_min > drag_limit
        # No panel angle can satisfy the structural limit. Command the
        # minimum-drag configuration as the protective fallback.
        return min_alpha
    end

    f(alpha) = q * _energy_depletion_struct_drag_area(
        spacecraft,
        temperature,
        speed_ratio,
        controlled_panel_links,
        alpha,
        config,
    ) - drag_limit

    alpha = try
        Roots.find_zero(f, (0.0, pi / 2), Roots.Bisection())
    catch
        min_alpha
    end
    return (isfinite(alpha) && 0.0 <= alpha <= max_alpha) ? alpha : 0.0
end

function _edg_structural_alpha(
    config::AerobrakingEnergyDepletionConfig,
    p::ODEParams,
    env,
    spacecraft,
    controlled_panel_links::Tuple{Vararg{Int}},
    base_alpha::Float64,
)
    _edg_in_drag_passage(p, env) || return base_alpha
    base_alpha < _edg_constraint_min_alpha(config) && return base_alpha
    return _energy_depletion_struct_load_root_alpha(config, env, spacecraft, controlled_panel_links, base_alpha)
end
