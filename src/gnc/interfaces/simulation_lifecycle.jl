# Internal lifecycle extension points. Legacy guidance retains no-op defaults.
module SimulationLifecycle
preflight_guidance(model, args; isolate_state) = nothing
initialize_guidance!(model, u, p, t) = nothing
before_reference_control!(guidance, control, u, p, t) = nothing
end
