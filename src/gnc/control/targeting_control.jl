using ..EDGAlgorithms
using ..EDGAlgorithms: _edg_orbit_metrics_from_rv,
    _edg_target_energy_from_apoapsis,
    _edg_targeting_constrained_alpha,
    _edg_targeting_prediction_time_grid,
    _edg_targeting_aero_acceleration,
    _edg_integrated_targeting_trajectory,
    _edg_integrated_max_energy_depletion_trajectory,
    _edg_predict_targeting_outcome,
    _edg_predict_max_energy_depletion_outcome,
    _edg_targeting_switch_outcomes,
    _edg_targeting_bracket_outcomes,
    _edg_targeting_outcome_with_heat_load,
    _edg_certify_targeting_candidates,
    _edg_disable_uncertified_targeting!,
    _edg_solve_targeting_switch
using ..EDGServices: _edg_control_sat_state,
    _edg_control_pos_vel_mass,
    _edg_environment_state,
    _edg_in_drag_passage,
    _edg_ephemeris_time,
    _edg_planet_frame_lpi,
    _edg_targeting_prediction_environment

"""
    SolarPanelAngleOfAttackControlModel(; controlled_panel_links=(2, 3))

Control effector that articulates the listed solar-panel links to realize a
commanded angle of attack during aerobraking passes. Link indices must be
positive and at least one link is required.
"""
struct SolarPanelAngleOfAttackControlModel <: AbstractControlEffectorModel
    controlled_panel_links::Tuple{Vararg{Int}}
end

@inline function _edg_panel_link_tuple(controlled_panel_links)
    links = tuple((Int(idx) for idx in controlled_panel_links)...)
    isempty(links) && throw(ArgumentError("SolarPanelAngleOfAttackControlModel requires at least one controlled link."))
    any(<=(0), links) && throw(ArgumentError("Controlled panel link indices must be positive."))
    return links
end

function SolarPanelAngleOfAttackControlModel(; controlled_panel_links=(2, 3))
    links = _edg_panel_link_tuple(controlled_panel_links)
    return SolarPanelAngleOfAttackControlModel(links)
end

"""
    AerobrakingEnergyDepletionControlModel

Control-side companion of the energy-depletion guidance strategy: tracks the
guidance-selected mode and drives the panel angle-of-attack effector under
the configured heat and structural limits. Holds the shared
[`AerobrakingEnergyDepletionConfig`](@ref) / [`AerobrakingEnergyDepletionState`](@ref)
and a [`SolarPanelAngleOfAttackControlModel`](@ref).
"""
struct AerobrakingEnergyDepletionControlModel <: AbstractControlEffectorModel
    config::AerobrakingEnergyDepletionConfig
    state::AerobrakingEnergyDepletionState
    aoa_effector::SolarPanelAngleOfAttackControlModel
end

function AerobrakingEnergyDepletionControlModel(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState;
    aoa_effector::SolarPanelAngleOfAttackControlModel=SolarPanelAngleOfAttackControlModel(config.controlled_panel_links),
)
    return AerobrakingEnergyDepletionControlModel(config, state, aoa_effector)
end

@inline function _edg_control_state_index_ok(state::AerobrakingEnergyDepletionState, i::Int)::Bool
    return 1 <= i <= length(state.selected_mode)
end

function _edg_recompute_switches!(
    model::AerobrakingEnergyDepletionControlModel,
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
    EDGAlgorithms._edg_recompute_switches!(model.config, model.state, p, env, spacecraft, pos, vel, mass, heat_load_j_cm2, t, i)
    return nothing
end

@inline function _edg_base_alpha(model::AerobrakingEnergyDepletionControlModel, t::Float64, i::Int)::Float64
    return EDGAlgorithms._edg_base_alpha(model.config, model.state, t, i)
end

function _edg_command_alpha!(
    model::AerobrakingEnergyDepletionControlModel,
    p::ODEParams,
    u,
    env,
    spacecraft,
    base_alpha::Float64,
    heat_load_j_cm2::Float64,
    heat_load_low_drag_active::Bool,
    i::Int,
)
    return EDGAlgorithms._edg_command_alpha!(model.config, model.state, p, model.aoa_effector.controlled_panel_links, env, spacecraft, base_alpha, heat_load_j_cm2, heat_load_low_drag_active, i)
end

function _apply_solar_panel_aoa!(
    effector::SolarPanelAngleOfAttackControlModel,
    spacecraft,
    alpha::Float64,
)
    links = spacecraft.links
    for idx in effector.controlled_panel_links
        1 <= idx <= length(links) || throw(ArgumentError("Controlled panel link index $(idx) is out of bounds for spacecraft with $(length(links)) links."))
        link = links[idx]
        link.root && throw(ArgumentError("Controlled panel link $(idx) is a root link and cannot be rotated."))
        axis = SVector{3, Float64}(abs.(link.r))
        rotate_link(link, axis, pi / 2 - alpha)
        link.α = alpha
    end
    return nothing
end

function calcControlEffect!(
    model::AerobrakingEnergyDepletionControlModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    decision = EDGAlgorithms.control_decision!(
        model.config, model.state, model.aoa_effector.controlled_panel_links, u, p, t, i,
    )
    decision.apply || return nothing
    _apply_solar_panel_aoa!(model.aoa_effector, p.args.dynamics_model.spacecraft[i], decision.command.alpha_command)
    return nothing
end

function calcControlForceTorque(
    model::AerobrakingEnergyDepletionControlModel,
    u::AbstractVector,
    p::ODEParams,
    i::Int64,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function calcControlForceTorque(
    model::SolarPanelAngleOfAttackControlModel,
    u::AbstractVector,
    p::ODEParams,
    i::Int64,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end
