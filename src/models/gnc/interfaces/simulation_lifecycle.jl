# Lifecycle hooks of guidance, navigation and control models: the pre-run checks,
# the run-start initializers, the atmosphere events and the call before a linked
# consumer's tick. Each default does nothing.
# The bind_spacecraft fallback names the effector's scope, a GNCRuntime declaration.
isdefined(@__MODULE__, :GNCRuntime) || Base.include(@__MODULE__, joinpath(@__DIR__, "gnc_runtime.jl"))
module SimulationLifecycle
using ..GNCRuntime: spacecraft_scope, BoundSpacecraft, CoupledSpacecraft

export preflight_guidance, preflight_control, initialize_navigation!, initialize_guidance!, initialize_control!,
    requires_atmosphere_events, atmosphere_transition!, before_reference_control!

"""
    preflight_guidance(model, args; isolate_state) -> nothing

Before a run starts, for each guidance `model` of the configuration `args`: check
that the configuration can run it and throw an `ArgumentError` naming the problem
otherwise. `isolate_state` is the run's setting. The check reads the caller's
objects, before the run copies them, so it must change nothing. Default: no check.
"""
preflight_guidance(model, args; isolate_state) = nothing

"""
    preflight_control(model, args; isolate_state) -> nothing

As `preflight_guidance`, for each control `model`. Default: no check.
"""
preflight_control(model, args; isolate_state) = nothing

"""
    initialize_navigation!(model, u, p, t) -> nothing

At the start of each run, for each navigation `model`: reset its run state, so
that two runs of one configuration with `isolate_state=false` start equal. `u` is
the initial state, `p` the run's parameters and `t` the start time (s into the
run). The engine calls the links' `link_initialize!` first, then this hook for
every navigation model, then `initialize_guidance!` and
`initialize_control!`, so each stage may read what the previous one set
up. What a run leaves in the model stays readable after the run: the reset
happens when the next run starts. Run state of a custom sensor is reset in its
`initialize_sensor!`. Default: nothing.
"""
initialize_navigation!(model, u, p, t) = nothing

"""
    initialize_guidance!(model, u, p, t) -> nothing

As `initialize_navigation!`, for each guidance `model`, after the navigation
models: reset its run state, so that two runs of one configuration with
`isolate_state=false` start equal. Every built-in guidance model's run state is
reset when a run starts, here or by the link it shares (`link_initialize!`).
Default: nothing.
"""
initialize_guidance!(model, u, p, t) = nothing

"""
    initialize_control!(model, u, p, t) -> nothing

As `initialize_navigation!`, for each control `model`, after the guidance
models: reset its run state, so that two runs of one configuration with
`isolate_state=false` start equal. Every built-in control model's run state is
reset when a run starts, here or by the link it shares; a value a user may
configure and the run then changes, such as a held command, returns to its value
before the model's first run. Default: nothing.
"""
initialize_control!(model, u, p, t) = nothing

"""
    requires_atmosphere_events(model) -> Bool

Whether guidance or control `model` needs `atmosphere_transition!` calls.
When any model of the run declares `true`, the engine tracks atmosphere entries
and exits even where its tolerances would not need them. Default: `false`.
"""
requires_atmosphere_events(model) = false

"""
    atmosphere_transition!(model, u, p, t, spacecraft_index, inside) -> nothing

Called for each guidance and control `model` when the spacecraft at
`spacecraft_index` (its position in the run, from 1) enters (`inside` is `true`)
or leaves the atmosphere at time `t` (s into the run), after the engine has
reconciled every spacecraft's crossing at that time. `u` is the state and `p` the
run's parameters. The model's run state is its own. Default: nothing.
"""
atmosphere_transition!(model, u, p, t, spacecraft_index, inside) = nothing

"""
    before_reference_control!(producer, consumer, u, p, t) -> nothing

Before each periodic tick of control model `consumer`, once per tick and on the
callback's thread, for each model that produces a link `consumer` consumes
(`gnc_links`), in configuration order: let `producer` bring what it writes up to
date for time `t` (s into the run), as the RPO planner checks its reference. `u`
is the state and `p` the run's parameters. Default: nothing.
"""
before_reference_control!(producer, consumer, u, p, t) = nothing

"""
    bind_spacecraft(effector, sat_idx::Int)

Extension point for per-spacecraft GNC: return `effector` bound to spacecraft
`sat_idx` (a run index, 1-based). A `SpacecraftModel` can declare guidance,
navigation and control effectors; at the start of `run_simulation` each one is
passed through this function with the position of its spacecraft, and the result
joins the configuration-level GNC tuples.

Implement it for an effector that selects its vehicle by an index field: return a
shallow copy (sharing mutable members with the original) with that field set, or
`effector` itself when it already targets `sat_idx`. A member that holds one
spacecraft's run state (a held command, a controller's warm start) needs its own
copy, so spacecraft bound from one template do not share it; members indexed by
spacecraft, such as an `ApolloDescentState`, stay shared. The fallback throws an
`ArgumentError` that names the effector's `spacecraft_scope`: an effector that acts
on every spacecraft (`AllSpacecraft`) or couples several vehicles
(`CoupledSpacecraft`, for example the RPO chaser and target types) cannot belong
to one; declare it at configuration level instead.
"""
function bind_spacecraft(effector, sat_idx::Int)
    T = typeof(effector)
    scope = spacecraft_scope(effector)
    why = scope isa CoupledSpacecraft ?
        "it couples several spacecraft (its spacecraft_scope is CoupledSpacecraft)" :
        scope isa BoundSpacecraft ?
        "it is bound to spacecraft $(join(scope.indices, ", ")) (its spacecraft_scope) and has no method to bind another" :
        "it acts on every spacecraft (its spacecraft_scope is AllSpacecraft)"
    throw(ArgumentError(
        "Cannot declare $(nameof(T)) on a spacecraft: $why. Declare it at configuration level " *
        "(guidance_model/navigation_model/control_model), or implement SpaceAGORA.bind_spacecraft(::$(nameof(T)), ::Int)."))
end

# Shallow copy of a plain struct with all fields set by position (kwdef types included).
_shallow_copy(x::T) where {T} = T((getfield(x, i) for i in 1:fieldcount(T))...)
end
