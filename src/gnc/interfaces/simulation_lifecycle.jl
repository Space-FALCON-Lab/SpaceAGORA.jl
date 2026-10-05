# Internal lifecycle extension points. Legacy guidance retains no-op defaults.
module SimulationLifecycle
preflight_guidance(model, args; isolate_state) = nothing
initialize_guidance!(model, u, p, t) = nothing
before_reference_control!(guidance, control, u, p, t) = nothing

"""
    bind_spacecraft(effector, sat_idx::Int)

Extension point for per-spacecraft GNC: return `effector` bound to spacecraft
`sat_idx` (a run index, 1-based). A `SpacecraftModel` can declare guidance,
navigation and control effectors; at the start of `run_simulation` each one is
passed through this function with the position of its spacecraft, and the result
joins the configuration-level GNC tuples.

Implement it for an effector that selects its vehicle by an index field: return a
shallow copy (sharing mutable members with the original) with that field set, or
`effector` itself when it already targets `sat_idx`. The fallback throws an
`ArgumentError`, because an effector without such a field acts on every
spacecraft and an effector that couples several vehicles (for example the RPO
chaser and target types) cannot belong to one; declare both at configuration
level instead.
"""
function bind_spacecraft(effector, sat_idx::Int)
    T = typeof(effector)
    why = hasfield(T, :chaser_idx) ?
        "it couples a chaser and a target spacecraft" :
        "it has no spacecraft index, so it acts on every spacecraft"
    throw(ArgumentError(
        "Cannot declare $(nameof(T)) on a spacecraft: $why. Declare it at configuration level " *
        "(guidance_model/navigation_model/control_model), or implement SpaceAGORA.bind_spacecraft(::$(nameof(T)), ::Int)."))
end

# Shallow copy of a plain struct with all fields set by position (kwdef types included).
_shallow_copy(x::T) where {T} = T((getfield(x, i) for i in 1:fieldcount(T))...)
end
