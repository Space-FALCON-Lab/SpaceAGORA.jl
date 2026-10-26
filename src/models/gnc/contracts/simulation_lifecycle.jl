# Lifecycle hooks of guidance, navigation and control models: the pre-run checks,
# the run-start initializers, the run state in checkpoint files, the atmosphere
# events and the call before a linked consumer's tick. Each default does nothing;
# restore_checkpoint_state! has none.
# The bind_spacecraft fallback names the effector's scope, a GNCRuntime declaration.
#
# The engine calls them in this order in a run (see "Run state: owners, resets
# and restores" in docs/contributing/architecture.md):
# 1. `bind_spacecraft` for each effector a spacecraft declares, then
#    `preflight_guidance` and `preflight_control`, on the caller's objects, before
#    the run copies them;
# 2. at the run's start, each link's `link_initialize!` (GNCRuntime), then
#    `initialize_navigation!`, `initialize_guidance!` and `initialize_control!`;
# 3. on a resumed run, `restore_checkpoint_state!`, after those initializers;
# 4. during the solve, the models' ticks, with `before_reference_control!` before
#    a linked consumer's periodic tick and `atmosphere_transition!` at atmosphere
#    crossings;
# 5. `checkpoint_state` at every checkpoint, and at the start of a checkpointed or
#    resumed run, to check that a model keeping run state can restore it.
isdefined(@__MODULE__, :GNCRuntime) || Base.include(@__MODULE__, joinpath(@__DIR__, "gnc_runtime.jl"))
module SimulationLifecycle
using ..GNCRuntime: spacecraft_scope, AllSpacecraft, BoundSpacecraft, CoupledSpacecraft

export preflight_guidance, preflight_control, initialize_navigation!, initialize_guidance!, initialize_control!,
    checkpoint_state, restore_checkpoint_state!,
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
reset when a run starts, here or by the link it shares (`link_initialize!`), and a
value a user may configure is kept when the user changed it, as for
`initialize_control!`. Default: nothing.
"""
initialize_guidance!(model, u, p, t) = nothing

"""
    initialize_control!(model, u, p, t) -> nothing

As `initialize_navigation!`, for each control `model`, after the guidance
models: reset its run state, so that two runs of one configuration with
`isolate_state=false` start equal. Every built-in control model's run state is
reset when a run starts, here or by the link it shares. A value a user may
configure that the run also writes, such as a held command, returns to its value
at the previous run's start while it still equals the model's last write in that
run; a value the user changed since is kept, and later runs start from it.
Default: nothing.
"""
initialize_control!(model, u, p, t) = nothing

"""
    checkpoint_state(model) -> data

The run state `model` keeps itself, for a checkpoint file: plain data (numbers,
arrays, named tuples, random streams), or `nothing` when it keeps none. The
engine calls it at every checkpoint for each guidance, navigation and control
model and each mounted sensor, and stores a copy of the data with the engine's
own run state; a resumed run hands it to `restore_checkpoint_state!`. Return
`nothing` for the whole run or data for the whole run. A model that returns data
implements `restore_checkpoint_state!`: a checkpointed run refuses at setup a
model that has the first and not the second. Default: `nothing`.
"""
checkpoint_state(model) = nothing

"""
    restore_checkpoint_state!(model, data) -> nothing

On a run resumed from a checkpoint file, put back the run state `data` that
`checkpoint_state(model)` returned when the file was written (read from the
file, so the model may keep it). For a guidance, navigation or control model
the engine calls it after the run-start initializers, so the restored state
replaces their reset; for a sensor before its sampling resumes (a resumed run
does not call `initialize_sensor!`). The resumed run refuses a checkpoint whose
models differ in kind or position from its own, or that holds no run state for a
model that keeps some. No default: a model whose `checkpoint_state` returns data
implements it.
"""
function restore_checkpoint_state! end

"""
    requires_atmosphere_events(model) -> Bool

Whether guidance or control `model` needs `atmosphere_transition!` calls.
When any model of the run declares `true`, the engine tracks atmosphere entries
and exits even where its tolerances would not need them. Default: `false`.
"""
requires_atmosphere_events(model) = false

"""
    atmosphere_transition!(model, u, p, t, spacecraft_index, inside) -> nothing

Called for each guidance and control `model` whose `spacecraft_scope` names the
spacecraft at `spacecraft_index` (its position in the run, from 1) when that
spacecraft enters (`inside` is `true`) or leaves the atmosphere at time `t` (s into
the run), after the engine has reconciled every spacecraft's crossing at that
time. SpaceAGORA 0.6 and earlier called every model. `u` is the state and `p` the
run's parameters. The model's run state is its own. Default: nothing.
"""
atmosphere_transition!(model, u, p, t, spacecraft_index, inside) = nothing

"""
    before_reference_control!(producer, consumer, u, p, t) -> nothing

Before each periodic tick of control model `consumer`, once per tick and on the
callback's thread, for each model that produces a link `consumer` consumes
(`gnc_links`): let `producer` bring what it writes up to date for time `t` (s into
the run), as the RPO planner checks its reference. `u` is the state and `p` the
run's parameters. Not called at a tick that runs `consumer` for no spacecraft, as
when every spacecraft in its scope is inactive and it declares `skip_inactive`.

The producers are called link by link, in the order the run's models first declare
those links, and within a link in configuration order. A producer of several of
`consumer`'s links is called once, for the first of them. This is not
configuration order: with A producing links X and Y, B producing Y and C producing
X, a consumer of X and Y gets A, C, B. A producer must not rely on another having
been called first. Default: nothing.
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
        scope isa AllSpacecraft ? "it acts on every spacecraft (its spacecraft_scope is AllSpacecraft)" :
        "it acts on the spacecraft its spacecraft_scope, a $(nameof(typeof(scope))), selects, and has no method to " *
        "bind one"
    # A model for every spacecraft keeps any per-spacecraft settings as one entry per
    # spacecraft, as the burn scheduler does.
    shared = !(scope isa AllSpacecraft) ? "" :
        ", as one model for every spacecraft with one entry per spacecraft in its per-spacecraft settings " *
        "(as the burn scheduler, BaseThrusterModel, is declared)"
    throw(ArgumentError(
        "Cannot declare $(nameof(T)) on a spacecraft: $why. Declare it at configuration level " *
        "(guidance_model/navigation_model/control_model)$shared, or implement " *
        "SpaceAGORA.bind_spacecraft(::$(nameof(T)), ::Int)."))
end

# Shallow copy of a plain struct with all fields set by position (kwdef types included).
_shallow_copy(x::T) where {T} = T((getfield(x, i) for i in 1:fieldcount(T))...)

# Values a user may configure that a run also writes, such as a held command or a
# plan buffer's plan. A model records its last write of them (`written`) during a
# run. At the next run's start, an entry still equal to that write returns to its
# value at the previous run's start (`start`, a copy); an entry changed since is the
# user's new configuration and is kept. The result is the run's start and the new
# baseline. Before the first run (`start === nothing`), or when the model wrote
# nothing, every entry is kept. Entries are the fields of `current`; `written` may
# hold more. The burn scheduler applies the same rule to its schedule entries.
function _run_start_entries(start, written, current::NamedTuple)
    (start === nothing || written === nothing) && return current
    names = keys(current)
    values = map(names) do name
        now = getfield(current, name)
        _same_value(now, getfield(written, name)) ? deepcopy(getfield(start, name)) : now
    end
    return NamedTuple{names}(values)
end

# Equal contents: numbers, strings, symbols and `nothing` by `isequal` (NaN equals
# NaN), arrays element by element, other objects field by field.
function _same_value(a, b)
    a === b && return true
    typeof(a) === typeof(b) || return false
    (a isa Number || a isa AbstractString || a isa Symbol) && return isequal(a, b)
    a isa AbstractArray && return size(a) == size(b) && all(i -> _same_value(a[i], b[i]), eachindex(a, b))
    a isa AbstractDict && return isequal(a, b)
    T = typeof(a)
    fieldcount(T) == 0 && return isequal(a, b)
    for i in 1:fieldcount(T)
        isdefined(a, i) == isdefined(b, i) || return false
        isdefined(a, i) && !_same_value(getfield(a, i), getfield(b, i)) && return false
    end
    return true
end
end
