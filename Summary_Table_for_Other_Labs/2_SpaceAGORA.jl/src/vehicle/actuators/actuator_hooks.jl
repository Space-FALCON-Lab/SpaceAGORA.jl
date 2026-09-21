module ActuatorHooks

export apply_actuator_hook!, apply_actuator_hooks!

"""
    apply_actuator_hook!(model, u, p, t, sat_idx)

Generic hook surface for actuator-side periodic updates.
Concrete actuator families should extend this method.
"""
function apply_actuator_hook!(model, u, p, t::Float64, sat_idx::Int)
    throw(MethodError(apply_actuator_hook!, (model, u, p, t, sat_idx)))
end

"""
    apply_actuator_hooks!(models, u, p, t, sat_idx)

Apply a tuple of actuator hooks for a single spacecraft index.
"""
function apply_actuator_hooks!(models::Tuple, u, p, t::Float64, sat_idx::Int)
    @inbounds for model in models
        apply_actuator_hook!(model, u, p, t, sat_idx)
    end
    return nothing
end

end # module ActuatorHooks
