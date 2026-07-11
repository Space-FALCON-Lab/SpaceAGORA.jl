using Roots

function _energy_depletion_struct_drag_area(
    spacecraft,
    temperature::Float64,
    speed_ratio::Float64,
    controlled_panel_links::Tuple{Vararg{Int}},
    controlled_alpha::Float64,
    config,
)::Float64
    controlled = Set{Int}(controlled_panel_links)
    drag_area = 0.0
    aero_module = getfield(getfield(parentmodule(@__MODULE__), :DynamicEffectors), :AerodynamicEffectors)
    aero_coeff = getfield(aero_module, :aerodynamic_coefficient_fM)

    for (idx, link) in pairs(spacecraft.links)
        area = max(0.0, Float64(link.ref_area))
        area == 0.0 && continue

        alpha = idx in controlled ?
            controlled_alpha :
            clamp(Float64(link.α), config.min_alpha_rad, config.max_alpha_rad)
        coeffs = aero_coeff(link, temperature, speed_ratio, alpha, Float64(link.β), Float64(link.θ))
        drag_area += max(0.0, Float64(coeffs[2])) * area
    end
    return drag_area
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

    min_alpha = config.min_alpha_rad
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
        return max_alpha
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
        Roots.find_zero(f, (min_alpha, max_alpha), Roots.Bisection())
    catch
        min_alpha
    end
    return clamp(alpha, min_alpha, max_alpha)
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
    return _energy_depletion_struct_load_root_alpha(config, env, spacecraft, controlled_panel_links, base_alpha)
end
