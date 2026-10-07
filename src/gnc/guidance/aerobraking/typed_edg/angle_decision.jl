function _edg_recompute_switches!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    p::ODEParams,
    env,
    spacecraft,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    mass::Float64,
    heat_load_j_cm2::Float64,
    t::Float64,
    i::Int,
)
    if state.selected_mode[i] == :targeting && state.targeting_active[i]
        isfinite(state.targeting_switch_s[i]) && return :cached
        _edg_in_drag_passage(p, env) || return :outside_pass
        state.targeting_switch_s[i] = _edg_solve_targeting_switch(
            config,
            state,
            p,
            spacecraft,
            pos,
            vel,
            mass,
            t,
            i;
            heat_load_j_cm2=heat_load_j_cm2,
            heat_rate_control=(:heat_rate in config.max_energy_submodes),
            structural_control=(:structural_load in config.max_energy_submodes),
        )
        state.last_switch_solve_t[i] = t
        return :solved
    elseif state.selected_mode[i] == :max_energy_depletion && (:heat_load in config.max_energy_submodes)
        if !_edg_in_drag_passage(p, env)
            if state.heat_load_drag_passage_active[i]
                state.heat_load_switch_solved[i] = false
                state.heat_load_drag_passage_active[i] = false
            end
            return :outside_pass
        end
        state.heat_load_drag_passage_active[i] = true
        state.heat_load_switch_solved[i] && return :cached
        state.heat_load_switches_s[i] = _edg_solve_heat_load_switches(
            config,
            p,
            spacecraft,
            pos,
            vel,
            mass,
            env,
            heat_load_j_cm2,
            t;
            heat_rate_control=(:heat_rate in config.max_energy_submodes),
            structural_control=(:structural_load in config.max_energy_submodes),
        )
        state.heat_load_switch_solved[i] = true
        state.last_switch_solve_t[i] = t
        return :solved
    end
    return :not_required
end

@inline function _edg_base_alpha(config::AerobrakingEnergyDepletionConfig, state::AerobrakingEnergyDepletionState, t::Float64, i::Int)::Float64
    mode = state.selected_mode[i]
    if mode == :safe_low_drag || state.safe_low_drag[i]
        return config.min_alpha_rad
    elseif mode == :targeting && state.targeting_active[i]
        return t >= state.targeting_switch_s[i] ? config.min_alpha_rad : config.max_alpha_rad
    elseif mode == :max_energy_depletion
        if _edg_heat_load_low_drag_active(config, state, t, i)
            return config.min_alpha_rad
        end
    end
    return config.max_alpha_rad
end

function _edg_command_alpha!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    p::ODEParams,
    controlled_panel_links::Tuple{Vararg{Int}},
    env,
    spacecraft,
    base_alpha::Float64,
    heat_load_j_cm2::Float64,
    heat_load_low_drag_active::Bool,
    i::Int,
)
    alpha = clamp(base_alpha, config.min_alpha_rad, config.max_alpha_rad)
    alpha_hr = alpha
    alpha_struct = alpha
    heat_rate_limit = config.heat_rate_limit_w_cm2
    heat_rate_active = (:heat_rate in config.max_energy_submodes) &&
        !heat_load_low_drag_active
    if heat_rate_active
        alpha_past = isfinite(state.last_alpha_rad[i]) ? state.last_alpha_rad[i] : alpha
        alpha_hr = _edg_heat_rate_alpha(config, p, env, alpha; limit_override=heat_rate_limit, alpha_past=alpha_past)
    end
    structural_active = (:structural_load in config.max_energy_submodes) &&
        !heat_load_low_drag_active
    if structural_active
        alpha_struct = _edg_structural_alpha(
            config,
            p,
            env,
            spacecraft,
            controlled_panel_links,
            alpha,
        )
    end
    if heat_rate_active &&
            structural_active
        alpha = min(alpha_hr, alpha_struct)
    elseif heat_rate_active
        alpha = alpha_hr
    elseif structural_active
        alpha = alpha_struct
    end
    if (:heat_load in config.max_energy_submodes) && heat_load_j_cm2 >= config.heat_load_limit_j_cm2
        alpha = min(alpha, config.min_alpha_rad)
    end
    alpha = clamp(alpha, config.min_alpha_rad, config.max_alpha_rad)
    state.last_alpha_rad[i] = alpha
    state.last_alpha_heat_rate_rad[i] = alpha_hr
    state.last_alpha_structural_rad[i] = alpha_struct
    state.last_heat_rate_w_cm2[i] = _edg_maxwellian_heat_rate(p, env, alpha)
    state.last_heat_load_j_cm2[i] = heat_load_j_cm2
    state.last_dynamic_pressure_pa[i] = env.dynamic_pressure
    return alpha
end

"""
Private result of a control-time EDG decision. The angle command is in radians;
apply distinguishes an invalid-index no-op from a valid zero angle. Diagnostics
reflect state already computed before actuation. switch_action describes cache/
pass/solve execution, not scientific acceptance of a solution.
"""
struct EDGControlDecision{D}
    command::AerobrakingControlCommand
    apply::Bool
    switch_action::Symbol
    diagnostics::D
end

function control_decision!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    controlled_panel_links::Tuple{Vararg{Int}},
    u,
    p::ODEParams,
    t::Float64,
    i::Int,
)
    1 <= i <= length(state.selected_mode) ||
        return EDGControlDecision(AerobrakingControlCommand(), false, :invalid_index, nothing)
    if state.selected_mode[i] == :inactive
        state.selected_mode[i] = (:max_energy_depletion in config.guidance_modes) ? :max_energy_depletion : :safe_low_drag
    end
    sc = _edg_control_sat_state(u, i)
    env = _edg_environment_state(u, p, Float64(t), i)
    spacecraft = p.args.dynamics_model.spacecraft[i]
    pos, vel, mass = _edg_control_pos_vel_mass(sc)
    heat_load = _edg_max_heat_load_for_links(sc, controlled_panel_links)
    switch_action = _edg_recompute_switches!(config, state, p, env, spacecraft, pos, vel, mass, heat_load, Float64(t), i)
    heat_load_low_drag_active = _edg_heat_load_low_drag_active(config, state, Float64(t), i)
    base_alpha = _edg_base_alpha(config, state, Float64(t), i)
    alpha = _edg_command_alpha!(config, state, p, controlled_panel_links, env, spacecraft, base_alpha, heat_load, heat_load_low_drag_active, i)
    diagnostics = (
        mode=state.selected_mode[i],
        alpha_heat_rate_rad=state.last_alpha_heat_rate_rad[i],
        alpha_structural_rad=state.last_alpha_structural_rad[i],
        heat_rate_w_cm2=state.last_heat_rate_w_cm2[i],
        heat_load_j_cm2=state.last_heat_load_j_cm2[i],
        dynamic_pressure_pa=state.last_dynamic_pressure_pa[i],
        targeting_switch_s=state.targeting_switch_s[i],
        heat_load_switches_s=state.heat_load_switches_s[i],
    )
    return EDGControlDecision(AerobrakingControlCommand(alpha), true, switch_action, diagnostics)
end
