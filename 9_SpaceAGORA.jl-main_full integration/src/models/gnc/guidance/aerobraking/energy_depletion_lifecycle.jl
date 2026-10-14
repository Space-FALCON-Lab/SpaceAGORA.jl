# The engine owns atmosphere crossings; EDG owns the lifetime of its cached plan.
function _edg_reset_pass!(s::AerobrakingEnergyDepletionState, i::Int)
    s.selected_mode[i] = :inactive
    s.targeting_active[i] = false
    s.safe_low_drag[i] = false
    s.energy_bracketing_evaluated[i] = false
    s.target_energy_jkg[i] = NaN
    s.bracket_min_energy_jkg[i] = NaN
    s.bracket_max_energy_jkg[i] = NaN
    s.heat_load_switches_s[i] = (Inf, Inf)
    s.heat_load_switch_solved[i] = false
    s.heat_load_drag_passage_active[i] = false
    s.targeting_switch_s[i] = Inf
    s.last_switch_solve_t[i] = -Inf
    s.heat_load_entry_time_s[i] = NaN
    s.heat_load_last_reevaluation_s[i] = -Inf
    s.heat_load_security_active[i] = false
    s.heat_load_previous_j_cm2[i] = 0.0
    s.target_perturbation_energy_change_jkg[i] = NaN
    s.targeting_diagnostics[i] = (status=:unplanned,)
    # Keep the cumulative solve count and last measured/commanded telemetry.
    # Propagated heat loads and actuator geometry are not EDG cache state.
    return nothing
end

function _edg_initialize_state!(s::AerobrakingEnergyDepletionState)
    s.control_clock_origin_s = 0.0
    for i in eachindex(s.selected_mode)
        _edg_reset_pass!(s, i)
        s.energy_bracketing_count[i] = 0
        s.last_alpha_rad[i] = NaN
        s.last_alpha_heat_rate_rad[i] = NaN
        s.last_alpha_structural_rad[i] = NaN
        s.last_heat_rate_w_cm2[i] = NaN
        s.last_heat_load_j_cm2[i] = NaN
        s.last_dynamic_pressure_pa[i] = NaN
        s.last_structural_load_pa[i] = NaN
        s.last_minimum_heat_rate_w_cm2[i] = NaN
        s.last_minimum_structural_load_pa[i] = NaN
        s.last_heat_rate_status[i] = :unobserved
        s.last_structural_load_status[i] = :unobserved
        s.last_heat_budget_status[i] = :unobserved
        s.heat_load_entry_j_cm2[i] = Float64[]
        s.heat_load_exit_j_cm2[i] = Float64[]
        s.last_pass_heat_load_j_cm2[i] = NaN
        s.heat_boundary_inside[i] = nothing
        s.heat_boundary_entry_time_s[i] = NaN
        s.heat_boundary_exit_time_s[i] = NaN
        s.heat_boundary_observed_at_s[i] = NaN
    end
    return nothing
end

function _edg_preflight(model, args)
    EDGServices.require_provider()
    n = length(args.dynamics_model.spacecraft)
    for field in fieldnames(AerobrakingEnergyDepletionState)
        field == :control_clock_origin_s && continue
        length(getfield(model.state, field)) == n || throw(ArgumentError(
            "EDG state.$field length must match the number of spacecraft ($n)."))
    end
    ss = args.simulation_settings
    (ss.checkpoint_enabled || ss.resume_from_checkpoint) && throw(ArgumentError(
        "EDG does not support checkpoint writing or resume: checkpoints do not contain its pass state."))
    guidance = filter(g -> g isa AerobrakingEnergyDepletionGuidanceModel,
        args.guidance_model.guidance_effectors)
    controls = filter(c -> c isa _control_module().AerobrakingEnergyDepletionControlModel,
        args.control_model.control_effectors)
    length(guidance) <= 1 && length(controls) <= 1 || throw(ArgumentError(
        "EDG requires at most one guidance model and one control model, with per-spacecraft state."))
    if !isempty(guidance) && !isempty(controls)
        g, c = only(guidance), only(controls)
        g.state === c.state || throw(ArgumentError("EDG guidance and control must share the same mutable state."))
        all(f -> isequal(getfield(g.config, f), getfield(c.config, f)),
            fieldnames(AerobrakingEnergyDepletionConfig)) || throw(ArgumentError(
            "EDG guidance and control configurations must agree."))
    elseif isempty(guidance) && :targeting in model.config.guidance_modes
        throw(ArgumentError("EDG targeting control requires its paired guidance model to bracket the target."))
    end
    return nothing
end

SimulationLifecycle.preflight_guidance(g::AerobrakingEnergyDepletionGuidanceModel, args; isolate_state) =
    _edg_preflight(g, args)
SimulationLifecycle.initialize_guidance!(g::AerobrakingEnergyDepletionGuidanceModel, u, p, t) =
    _edg_initialize_state!(g.state)
SimulationLifecycle.atmosphere_transition!(g::AerobrakingEnergyDepletionGuidanceModel, u, p, t, i, inside) =
    nothing # EDG resets at its own geodetic passage surface.

SimulationLifecycle.requires_atmosphere_events(::AerobrakingEnergyDepletionGuidanceModel) = true
