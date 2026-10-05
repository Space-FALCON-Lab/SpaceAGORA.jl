# Internal lifecycle extension points. Legacy guidance retains no-op defaults.
module SimulationLifecycle
preflight_guidance(model, args; isolate_state) = nothing
initialize_guidance!(model, u, p, t) = nothing
preflight_control(model, args; isolate_state) = nothing
initialize_control!(model, u, p, t) = nothing
# Called after the engine reconciles the full atmosphere mask, once per changed
# spacecraft. Model methods own their state; legacy models remain inert.
requires_atmosphere_events(model) = false
atmosphere_transition!(model, u, p, t, i, inside) = nothing
before_reference_control!(guidance, control, u, p, t) = nothing
end
