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
    # Keep the cumulative solve count and last measured/commanded telemetry.
    # Propagated heat loads and actuator geometry are not EDG cache state.
    return nothing
end

function _edg_initialize_state!(s::AerobrakingEnergyDepletionState)
    for i in eachindex(s.selected_mode)
        _edg_reset_pass!(s, i)
        s.energy_bracketing_count[i] = 0
        s.last_alpha_rad[i] = NaN
        s.last_alpha_heat_rate_rad[i] = NaN
        s.last_alpha_structural_rad[i] = NaN
        s.last_heat_rate_w_cm2[i] = NaN
        s.last_heat_load_j_cm2[i] = NaN
        s.last_dynamic_pressure_pa[i] = NaN
    end
    return nothing
end

function _edg_preflight(model, args)
    n = length(args.dynamics_model.spacecraft)
    for field in fieldnames(AerobrakingEnergyDepletionState)
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
    _edg_reset_pass!(g.state, i)

SimulationLifecycle.requires_atmosphere_events(::AerobrakingEnergyDepletionGuidanceModel) = true
