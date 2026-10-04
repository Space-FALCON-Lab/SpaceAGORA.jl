import ..GuidanceHooks

# Paired guidance owns initialization and crossing resets. A maximum-depletion
# controller also works alone, and then owns those operations itself.
_edg_has_guidance_owner(c, args) = any(g ->
    g isa GuidanceHooks.AerobrakingEnergyDepletionGuidanceModel && g.state === c.state,
    args.guidance_model.guidance_effectors)

SimulationLifecycle.preflight_control(c::AerobrakingEnergyDepletionControlModel, args; isolate_state) =
    GuidanceHooks._edg_preflight(c, args)
function SimulationLifecycle.initialize_control!(c::AerobrakingEnergyDepletionControlModel, u, p, t)
    _edg_has_guidance_owner(c, p.args) || GuidanceHooks._edg_initialize_state!(c.state)
    return nothing
end
function SimulationLifecycle.atmosphere_transition!(c::AerobrakingEnergyDepletionControlModel, u, p, t, i, inside)
    _edg_has_guidance_owner(c, p.args) || GuidanceHooks._edg_reset_pass!(c.state, i)
    return nothing
end

SimulationLifecycle.requires_atmosphere_events(::AerobrakingEnergyDepletionControlModel) = true
