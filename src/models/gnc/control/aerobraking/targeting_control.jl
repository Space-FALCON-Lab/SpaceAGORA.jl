using ..EDGAlgorithms
using ..EDGServices: _edg_control_sat_state,
    _edg_control_pos_vel_mass,
    _edg_ephemeris_time

"""
    SolarPanelAngleOfAttackControlModel(; controlled_panel_links=(2, 3))

Control effector that articulates the listed solar-panel links to realize a
commanded angle of attack during aerobraking passes. Link indices must be
positive and at least one link is required.

A panel link turns at once to the commanded angle of attack α (a rotation by
π/2 - α about the absolute value of its offset `r`), unless a `PanelHinge` is
mounted on it with `add_actuator!`: then the hinge moves toward π/2 - α at
each control tick within its angle and rate limits, and the panel's angle of
attack is the one it reached. With the hinge's axis along the absolute value
of the link's offset and no limits, the panel moves exactly as without one.
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

function _apply_solar_panel_aoa!(
    effector::SolarPanelAngleOfAttackControlModel,
    spacecraft,
    alpha::Float64,
    t::Float64=NaN,
)
    links = spacecraft.links
    for idx in effector.controlled_panel_links
        1 <= idx <= length(links) || throw(ArgumentError("Controlled panel link index $(idx) is out of bounds for spacecraft with $(length(links)) links."))
        link = links[idx]
        link.root && throw(ArgumentError("Controlled panel link $(idx) is a root link and cannot be rotated."))
        if isempty(link.hinges)
            axis = SVector{3, Float64}(abs.(link.r))
            rotate_link(link, axis, pi / 2 - alpha)
            link.α = alpha
        else
            # A mounted hinge moves the panel within its limits; the angle of
            # attack is the one it reached (exactly `alpha` when it got there).
            isfinite(t) || throw(ArgumentError(
                "Controlled panel link $(idx) carries a hinge, which moves at control ticks: pass the tick time."))
            target = pi / 2 - alpha
            angle = move_hinge!(link, link.hinges[1], target, t)
            link.α = angle == target ? alpha : pi / 2 - angle
        end
    end
    return nothing
end

function update_control!(
    model::AerobrakingEnergyDepletionControlModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    decision = EDGAlgorithms.control_decision!(
        model.state, model.config, model.aoa_effector.controlled_panel_links, u, p, t, i,
    )
    decision.apply || return nothing
    _apply_solar_panel_aoa!(model.aoa_effector, p.args.dynamics_model.spacecraft[i], decision.command.alpha_command, t)
    return nothing
end

function control_wrench(
    model::AerobrakingEnergyDepletionControlModel,
    x::StateSample,
    p::ODEParams,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

function control_wrench(
    model::SolarPanelAngleOfAttackControlModel,
    x::StateSample,
    p::ODEParams,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end


@inline function _edg_control_state_index_ok(state::AerobrakingEnergyDepletionState, i::Int)::Bool
    return 1 <= i <= length(state.selected_mode)
end

"""
    apply_panel_angle_of_attack!(effector::SolarPanelAngleOfAttackControlModel, spacecraft, alpha_rad, t)

Apply a panel angle of attack in radians at control-tick time `t` (seconds).
Only the effector's selected non-root links are moved. Hinge-less panels rotate
directly; mounted hinges respect their configured angle and rate limits, and
`link.α` records the reached angle. Invalid link indices and root links throw.
A mounted hinge requires a finite tick time. Call from committed control ticks,
not from a force evaluation or an optimizer's trial state.
"""
apply_panel_angle_of_attack!(effector::SolarPanelAngleOfAttackControlModel,
    spacecraft, alpha_rad::Real, t::Real) =
    _apply_solar_panel_aoa!(effector, spacecraft, Float64(alpha_rad), Float64(t))
