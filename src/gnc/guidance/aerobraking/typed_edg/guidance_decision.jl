function _edg_interpolate_bracket_value(exit_energy::Float64, energy_min::Float64, energy_max::Float64, value_at_min::Float64, value_at_max::Float64)
    width = energy_max - energy_min
    abs(width) < eps(Float64) && return 0.5 * (value_at_min + value_at_max)
    fraction = (exit_energy - energy_min) / width
    return value_at_min + fraction * (value_at_max - value_at_min)
end

function _edg_target_energy_from_reachable_bracket(
    planet,
    target_apoapsis_radius_m::Float64,
    energy_min::Float64,
    energy_max::Float64,
    periapsis_at_min::Float64,
    periapsis_at_max::Float64,
)
    function residual(exit_energy)
        periapsis = _edg_interpolate_bracket_value(exit_energy, energy_min, energy_max, periapsis_at_min, periapsis_at_max)
        desired_energy = _edg_target_energy_from_apoapsis(planet, target_apoapsis_radius_m, periapsis)
        return exit_energy - desired_energy
    end

    residual_min = residual(energy_min)
    residual_max = residual(energy_max)
    if isfinite(residual_min) && isfinite(residual_max) && residual_min * residual_max <= 0.0
        return Roots.find_zero(residual, (energy_min, energy_max), Roots.Brent(); rtol=1e-10)
    elseif isfinite(residual_min) && isfinite(residual_max) && abs(residual_max - residual_min) > eps(Float64)
        return energy_min - residual_min * (energy_max - energy_min) / (residual_max - residual_min)
    end
    return abs(residual_min) <= abs(residual_max) ? energy_min : energy_max
end

function _edg_set_targeting_fallback!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    i::Int,
)
    state.targeting_active[i] = false
    state.safe_low_drag[i] = !(:max_energy_depletion in config.guidance_modes)
    state.selected_mode[i] = (:max_energy_depletion in config.guidance_modes) ? :max_energy_depletion : :safe_low_drag
    return nothing
end

function _edg_run_target_energy_bracketing!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    u,
    p::ODEParams,
    t::Float64,
    i::Int,
)
    if state.energy_bracketing_evaluated[i]
        return nothing
    end

    env = _edg_environment_state(u, p, t, i)
    if !_edg_in_drag_passage(p, env)
        _edg_set_targeting_fallback!(config, state, i)
        return nothing
    end

    sc = _edg_control_sat_state(u, i)
    pos, vel, mass = _edg_control_pos_vel_mass(sc)
    spacecraft = p.args.dynamics_model.spacecraft[i]
    planet = p.args.environment_model.planet
    heat_load = _edg_max_heat_load_for_links(sc, config.controlled_panel_links)

    low_drag, max_energy_depletion = _edg_targeting_bracket_outcomes(
        config,
        p,
        spacecraft,
        pos,
        vel,
        mass,
        t;
        heat_load_j_cm2=heat_load,
        heat_rate_control=(:heat_rate in config.max_energy_submodes),
        structural_control=(:structural_load in config.max_energy_submodes),
    )

    endpoints = (low_drag, max_energy_depletion)
    energy_values = (low_drag.energy_jkg, max_energy_depletion.energy_jkg)
    energy_min, energy_max = extrema(energy_values)
    min_idx = energy_values[1] <= energy_values[2] ? 1 : 2
    max_idx = min_idx == 1 ? 2 : 1
    periapsis_at_min = endpoints[min_idx].periapsis_radius_m
    periapsis_at_max = endpoints[max_idx].periapsis_radius_m
    apoapsis_min, apoapsis_max = extrema((low_drag.apoapsis_radius_m, max_energy_depletion.apoapsis_radius_m))

    target_energy = _edg_target_energy_from_reachable_bracket(
        planet,
        config.target_apoapsis_radius_m,
        energy_min,
        energy_max,
        periapsis_at_min,
        periapsis_at_max,
    )

    energy_tol = 1e-6 * max(abs(energy_min), abs(energy_max), 1.0)
    apo_tol = 1e-6 * max(abs(apoapsis_min), abs(apoapsis_max), abs(config.target_apoapsis_radius_m), 1.0)
    reachable = isfinite(target_energy) &&
        energy_min - energy_tol <= target_energy <= energy_max + energy_tol &&
        apoapsis_min - apo_tol <= config.target_apoapsis_radius_m <= apoapsis_max + apo_tol

    state.energy_bracketing_evaluated[i] = true
    state.energy_bracketing_count[i] += 1
    state.target_energy_jkg[i] = target_energy
    state.bracket_min_energy_jkg[i] = energy_min
    state.bracket_max_energy_jkg[i] = energy_max
    state.targeting_active[i] = reachable
    state.safe_low_drag[i] = !reachable && !(:max_energy_depletion in config.guidance_modes)
    state.selected_mode[i] = reachable ? :targeting :
        ((:max_energy_depletion in config.guidance_modes) ? :max_energy_depletion : :safe_low_drag)
    return nothing
end

function guidance_decision!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    (1 <= i <= length(state.selected_mode)) || return nothing

    if :targeting in config.guidance_modes
        _edg_run_target_energy_bracketing!(config, state, u, p, Float64(t), i)
    elseif :max_energy_depletion in config.guidance_modes
        state.selected_mode[i] = :max_energy_depletion
        state.targeting_active[i] = false
        state.safe_low_drag[i] = false
        state.energy_bracketing_evaluated[i] = false
    else
        state.selected_mode[i] = :safe_low_drag
        state.targeting_active[i] = false
        state.safe_low_drag[i] = true
    end
    return nothing
end
