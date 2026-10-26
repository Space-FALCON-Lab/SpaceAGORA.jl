# Campaign supervisor: maximum energy depletion until the target energy is
# reachable in the coming pass, then energy targeting.
"""
Run state of an [`AerobrakingMPCCampaignControlModel`](@ref): its `phase`, the
passes completed (atmosphere exits of the run), the latest reachable-energy
bracket (MJ/kg) and the switch time.
"""
Base.@kwdef mutable struct AerobrakingMPCCampaignState
    phase::Symbol = :maximum_energy_depletion
    completed_passes::Int = 0
    bracket_min_energy_mj_kg::Float64 = NaN
    bracket_max_energy_mj_kg::Float64 = NaN
    switch_time_s::Float64 = NaN
    last_error::Union{Nothing, String} = nothing
end

"""
Aerobraking campaign supervisor: runs `maximum_depletion_control` until at
least `minimum_depletion_passes` passes are complete and the target energy lies
within the energy reachable in the coming pass, then hands over to
`targeting_control`. Both controllers must act on the same spacecraft and panels;
`run_simulation` checks the pair before the run starts.
"""
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

# Two-body specific energy (MJ/kg) of a KS state.
function _mpc_specific_energy_mj_kg(state, params::KSPropagationParams)
    cart = ks_state_to_cartesian(state)
    return (0.5 * dot(cart.velocity_ii_m, cart.velocity_ii_m) -
        params.μ / norm(cart.position_ii_m)) / 1.0e6
end

# The end-of-pass energy reachable from the current state lies between the
# maximum-depletion prediction (fully exposed panels, as planned) and a reference
# pass with the panels stowed (bus area only). Records that bracket and returns
# whether the target energy lies within it, up to `bracket_tolerance_mj_kg`.
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
        epoch_s=t,
    )
    rotation = _mpc_planet_frame_rotation(p, t)
    low_area = med.config.bus_reference_area_m2
    low_reference = build_reference_drag_pass(
        params,
        rotation * pos,
        rotation * vel;
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
    target_energy = model.targeting_control.config.target_energy_mj_kg
    tol = model.bracket_tolerance_mj_kg
    return model.state.bracket_min_energy_mj_kg - tol <= target_energy <=
        model.state.bracket_max_energy_mj_kg + tol
end

# Hand the held command and this pass's heat accounting to the targeting
# controller, which then solves at once.
function _mpc_campaign_transfer_state!(target::AerobrakingMPCState, source::AerobrakingMPCState)
    target.held_commanded_area_m2 = source.held_commanded_area_m2
    target.held_alpha_rad = source.held_alpha_rad
    target.last_command_time_s = source.last_command_time_s
    target.estimated_heat_load_j_cm2 = source.estimated_heat_load_j_cm2
    target.last_heat_rate_w_cm2 = source.last_heat_rate_w_cm2
    target.last_drag_n = source.last_drag_n
    target.last_solve_time_s = -Inf
    return nothing
end

function update_control!(
    model::AerobrakingMPCCampaignControlModel,
    u,
    p,
    t::Float64,
    sat_idx::Int64,
)
    med = model.maximum_depletion_control
    Int(sat_idx) == med.spacecraft_index || return nothing
    if model.state.phase == :maximum_energy_depletion
        previous_solve_count = med.state.solve_count
        update_control!(med, u, p, t, sat_idx)
        solved_now = med.state.solve_count > previous_solve_count
        if solved_now && model.state.completed_passes >= model.minimum_depletion_passes
            try
                if _mpc_campaign_reachable_bracket!(model, u, p, t, Int(sat_idx))
                    _mpc_campaign_transfer_state!(model.targeting_control.state, med.state)
                    model.state.phase = :target_energy
                    model.state.switch_time_s = t
                    model.state.last_error = nothing
                    update_control!(model.targeting_control, u, p, t, sat_idx)
                end
            catch err
                model.state.last_error = sprint(showerror, err)
            end
        end
    else
        update_control!(model.targeting_control, u, p, t, sat_idx)
    end
    return nothing
end

function control_wrench(
    model::AerobrakingMPCCampaignControlModel,
    x::StateSample,
    p,
    t::Float64,
)
    return SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)
end

_mpc_active_control(model::AerobrakingMPCCampaignControlModel) =
    model.state.phase == :target_energy ? model.targeting_control : model.maximum_depletion_control

"""
    mpc_campaign_save_fields(model) -> Tuple{SaveField...}

The [`mpc_control_save_fields`](@ref) columns of each campaign's active
controller plus its phase code (1 maximum depletion, 2 energy targeting),
completed passes and reachable-energy bracket (MJ/kg), from the run's own
[`AerobrakingMPCCampaignControlModel`](@ref)s.
"""
function mpc_campaign_save_fields(::AerobrakingMPCCampaignControlModel)
    T = AerobrakingMPCCampaignControlModel
    return (
        _mpc_control_columns(T)...,
        _mpc_save_field(:mpc_campaign_phase_code, T, e -> e.state.phase == :target_energy ? 2.0 : 1.0),
        _mpc_save_field(:mpc_campaign_completed_passes, T, e -> e.state.completed_passes),
        _mpc_save_field(:mpc_campaign_bracket_min_energy_mj_kg, T, e -> e.state.bracket_min_energy_mj_kg),
        _mpc_save_field(:mpc_campaign_bracket_max_energy_mj_kg, T, e -> e.state.bracket_max_energy_mj_kg),
    )
end

GNCRuntime.spacecraft_scope(m::AerobrakingMPCCampaignControlModel) =
    BoundSpacecraft((m.maximum_depletion_control.spacecraft_index,))
GNCRuntime.held_command(m::AerobrakingMPCCampaignControlModel, spacecraft_index) =
    GNCRuntime.held_command(_mpc_active_control(m), spacecraft_index)
GNCRuntime.held_command_time(m::AerobrakingMPCCampaignControlModel, spacecraft_index) =
    GNCRuntime.held_command_time(_mpc_active_control(m), spacecraft_index)

function SimulationLifecycle.bind_spacecraft(m::AerobrakingMPCCampaignControlModel, sat_idx::Int)
    m.maximum_depletion_control.spacecraft_index == sat_idx && return m
    b = SimulationLifecycle._shallow_copy(m)
    b.maximum_depletion_control = SimulationLifecycle.bind_spacecraft(m.maximum_depletion_control, sat_idx)
    b.targeting_control = SimulationLifecycle.bind_spacecraft(m.targeting_control, sat_idx)
    b.state = deepcopy(m.state)
    return b
end

SimulationLifecycle.requires_atmosphere_events(::AerobrakingMPCCampaignControlModel) = true
function SimulationLifecycle.atmosphere_transition!(
    m::AerobrakingMPCCampaignControlModel, u, p, t, spacecraft_index, inside)
    # A pass is complete when the spacecraft leaves the atmosphere.
    !inside && spacecraft_index == m.maximum_depletion_control.spacecraft_index &&
        (m.state.completed_passes += 1)
    SimulationLifecycle.atmosphere_transition!(m.maximum_depletion_control, u, p, t, spacecraft_index, inside)
    SimulationLifecycle.atmosphere_transition!(m.targeting_control, u, p, t, spacecraft_index, inside)
    return nothing
end

SimulationLifecycle.preflight_control(m::AerobrakingMPCCampaignControlModel, args; isolate_state) =
    _validate_mpc_campaign(m)
