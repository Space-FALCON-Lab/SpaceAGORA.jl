#=
"""
    Multi-pass aerobraking MPC campaign supervisor.

    The supervisor runs maximum-energy-depletion MPC for complete atmospheric
    passes. Before each subsequent pass it predicts the reachable terminal
    energy interval between minimum exposed area and the constrained MED plan.
    Once the user target lies in that interval, it switches permanently to the
    terminal-energy MPC controller.
"""
=#
Base.@kwdef mutable struct AerobrakingMPCCampaignState
    phase::Symbol = :maximum_energy_depletion
    completed_passes::Int = 0
    inside_atmosphere::Bool = false
    bracket_min_energy_mj_kg::Float64 = NaN
    bracket_max_energy_mj_kg::Float64 = NaN
    bracket_evaluations::Int = 0
    switch_time_s::Float64 = NaN
    last_bracket_time_s::Float64 = -Inf
    last_error::Union{Nothing, String} = nothing
end

Base.@kwdef mutable struct AerobrakingMPCCampaignControlModel <: AbstractControlEffectorModel
    maximum_depletion_control::AerobrakingMPCControlModel
    targeting_control::AerobrakingMPCControlModel
    state::AerobrakingMPCCampaignState = AerobrakingMPCCampaignState()
    minimum_depletion_passes::Int = 2
    bracket_tolerance_mj_kg::Float64 = 0.0
end

function _validate_mpc_campaign(model::AerobrakingMPCCampaignControlModel)
    med = model.maximum_depletion_control
    target = model.targeting_control
    med.config.mode isa MaxEnergyDepletionMode ||
        throw(ArgumentError("maximum_depletion_control must use MaxEnergyDepletionMode()."))
    target.config.mode isa TargetEnergyMode ||
        throw(ArgumentError("targeting_control must use TargetEnergyMode()."))
    med.spacecraft_index == target.spacecraft_index ||
        throw(ArgumentError("Campaign MPC controllers must use the same spacecraft_index."))
    med.controlled_panel_links == target.controlled_panel_links ||
        throw(ArgumentError("Campaign MPC controllers must use the same controlled_panel_links."))
    model.minimum_depletion_passes >= 1 ||
        throw(ArgumentError("minimum_depletion_passes must be at least one."))
    model.bracket_tolerance_mj_kg >= 0.0 ||
        throw(ArgumentError("bracket_tolerance_mj_kg must be nonnegative."))
    return nothing
end

function _mpc_specific_energy_mj_kg(state, params::AerobrakingMPCParams)
    cart = ks_state_to_cartesian(state)
    return (0.5 * dot(cart.velocity_ii_m, cart.velocity_ii_m) -
        params.μ / norm(cart.position_ii_m)) / 1.0e6
end

function _mpc_campaign_update_pass_count!(
    model::AerobrakingMPCCampaignControlModel,
    pos::SVector{3, Float64},
    vel::SVector{3, Float64},
    params::AerobrakingMPCParams,
)
    h = norm(pos) - params.Re
    radial_velocity = dot(pos, vel) / norm(pos)
    inside = h <= model.maximum_depletion_control.reference.h_cut_m
    if model.state.inside_atmosphere && !inside && radial_velocity > 0.0
        model.state.completed_passes += 1
    end
    model.state.inside_atmosphere = inside
    return nothing
end

function _mpc_campaign_reachable_bracket!(
    model::AerobrakingMPCCampaignControlModel,
    u,
    p,
    t::Float64,
    sat_idx::Int,
)
    med = model.maximum_depletion_control
    sc_state = _mpc_control_sat_state(u, sat_idx)
    pos, vel = _mpc_control_pos_vel(sc_state)
    params = mpc_params_from_spaceagora(p.args)
    density = density_function_from_spaceagora(
        p;
        latitude=med.prediction_latitude_rad,
        longitude=med.prediction_longitude_rad,
        wind=med.prediction_wind,
    )
    low_area = med.config.bus_reference_area_m2
    low_reference = build_reference_drag_pass(
        params,
        pos,
        vel;
        config=med.config,
        reference=med.reference,
        nominal_area_m2=low_area,
        density=density,
    )
    low_energy = _mpc_specific_energy_mj_kg(view(low_reference.states, size(low_reference.states, 1), :), params)
    med_energy = med.state.predicted_terminal_energy / 1.0e6
    isfinite(med_energy) || return false

    model.state.bracket_min_energy_mj_kg = min(low_energy, med_energy)
    model.state.bracket_max_energy_mj_kg = max(low_energy, med_energy)
    model.state.bracket_evaluations += 1
    model.state.last_bracket_time_s = t
    target_energy = model.targeting_control.config.target_energy_mj_kg
    tol = model.bracket_tolerance_mj_kg
    return model.state.bracket_min_energy_mj_kg - tol <= target_energy <=
        model.state.bracket_max_energy_mj_kg + tol
end

function _mpc_campaign_transfer_state!(target::AerobrakingMPCState, source::AerobrakingMPCState)
    target.held_commanded_area_m2 = source.held_commanded_area_m2
    target.held_alpha_rad = source.held_alpha_rad
    target.last_command_time_s = source.last_command_time_s
    target.estimated_heat_load_j_cm2 = source.estimated_heat_load_j_cm2
    target.last_heat_rate_w_cm2 = source.last_heat_rate_w_cm2
    target.last_drag_n = source.last_drag_n
    target.last_density_kg_m3 = source.last_density_kg_m3
    target.last_solve_time_s = -Inf
    return nothing
end

function calcControlEffect!(
    model::AerobrakingMPCCampaignControlModel,
    u,
    p,
    t::Float64,
    sat_idx::Int64,
)
    _validate_mpc_campaign(model)
    med = model.maximum_depletion_control
    Int(sat_idx) == med.spacecraft_index || return nothing
    sc_state = _mpc_control_sat_state(u, Int(sat_idx))
    pos, vel = _mpc_control_pos_vel(sc_state)
    params = mpc_params_from_spaceagora(p.args)
    _mpc_campaign_update_pass_count!(model, pos, vel, params)

    if model.state.phase == :maximum_energy_depletion
        previous_solve_count = med.state.solve_count
        calcControlEffect!(med, u, p, t, sat_idx)
        solved_now = med.state.solve_count > previous_solve_count
        if solved_now && model.state.completed_passes >= model.minimum_depletion_passes
            try
                if _mpc_campaign_reachable_bracket!(model, u, p, t, Int(sat_idx))
                    _mpc_campaign_transfer_state!(model.targeting_control.state, med.state)
                    model.state.phase = :target_energy
                    model.state.switch_time_s = t
                    model.state.last_error = nothing
                    calcControlEffect!(model.targeting_control, u, p, t, sat_idx)
                end
            catch err
                model.state.last_error = sprint(showerror, err)
            end
        end
    else
        calcControlEffect!(model.targeting_control, u, p, t, sat_idx)
    end
    return nothing
end

function calcControlForceTorque(
    model::AerobrakingMPCCampaignControlModel,
    u::AbstractVector,
    p,
    i::Int64,
    t::Float64,
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

_mpc_campaign_active_control(model::AerobrakingMPCCampaignControlModel) =
    model.state.phase == :target_energy ? model.targeting_control : model.maximum_depletion_control

function mpc_campaign_save_fields(model::AerobrakingMPCCampaignControlModel)
    SaveField = getproperty(parentmodule(@__MODULE__), :SaveField)
    satellite_count = integrator -> eachindex(integrator.p.args.dynamics_model.spacecraft)
    per_satellite_value = (integrator, value) -> [
        i == model.maximum_depletion_control.spacecraft_index ? value : NaN
        for i in satellite_count(integrator)
    ]
    return (
        SaveField(:mpc_commanded_area_m2, (u, t, integrator) ->
            per_satellite_value(integrator, _mpc_save_control(_mpc_campaign_active_control(model))[1]);
            per_satellite=true, column_prefix="mpc_commanded_area_m2"),
        SaveField(:mpc_alpha_rad, (u, t, integrator) ->
            per_satellite_value(integrator, _mpc_save_control(_mpc_campaign_active_control(model))[2]);
            per_satellite=true, column_prefix="mpc_alpha_rad"),
        SaveField(:mpc_active_limit_code, (u, t, integrator) ->
            per_satellite_value(integrator, _mpc_limit_code(_mpc_campaign_active_control(model).state.active_limit));
            per_satellite=true, column_prefix="mpc_active_limit_code"),
        SaveField(:mpc_heat_rate_w_cm2, (u, t, integrator) ->
            per_satellite_value(integrator, _mpc_campaign_active_control(model).state.last_heat_rate_w_cm2);
            per_satellite=true, column_prefix="mpc_heat_rate_w_cm2"),
        SaveField(:mpc_drag_n, (u, t, integrator) ->
            per_satellite_value(integrator, _mpc_campaign_active_control(model).state.last_drag_n);
            per_satellite=true, column_prefix="mpc_drag_n"),
        SaveField(:mpc_heat_load_j_cm2, (u, t, integrator) ->
            per_satellite_value(integrator, _mpc_campaign_active_control(model).state.estimated_heat_load_j_cm2);
            per_satellite=true, column_prefix="mpc_heat_load_j_cm2"),
        SaveField(:mpc_campaign_phase_code, (u, t, integrator) ->
            per_satellite_value(integrator, model.state.phase == :target_energy ? 2.0 : 1.0);
            per_satellite=true, column_prefix="mpc_campaign_phase_code"),
        SaveField(:mpc_campaign_completed_passes, (u, t, integrator) ->
            per_satellite_value(integrator, Float64(model.state.completed_passes));
            per_satellite=true, column_prefix="mpc_campaign_completed_passes"),
        SaveField(:mpc_campaign_bracket_min_energy_mj_kg, (u, t, integrator) ->
            per_satellite_value(integrator, model.state.bracket_min_energy_mj_kg);
            per_satellite=true, column_prefix="mpc_campaign_bracket_min_energy_mj_kg"),
        SaveField(:mpc_campaign_bracket_max_energy_mj_kg, (u, t, integrator) ->
            per_satellite_value(integrator, model.state.bracket_max_energy_mj_kg);
            per_satellite=true, column_prefix="mpc_campaign_bracket_max_energy_mj_kg"),
    )
end
