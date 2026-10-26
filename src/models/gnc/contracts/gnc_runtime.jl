# Declarations a control model makes about how the engine runs it. Each has a
# default that keeps the engine's usual behaviour, so a model declares only what
# differs from a periodic tick.
module GNCRuntime

export tick_schedule, gnc_callbacks, register_tstops!, gnc_clock_started!, gnc_period_s
export Produces, Consumes, gnc_links, link_cardinality, check_link, link_initialize!
export SpacecraftScope, AllSpacecraft, BoundSpacecraft, CoupledSpacecraft, spacecraft_scope, skip_inactive,
    held_command, held_command_time, release_command!

"""
    tick_schedule(model) -> Symbol

How the engine runs control `model`.

- `:periodic`, the default: a periodic callback at the model's control rate calls
  its tick (`update_control!`) for every spacecraft.
- `:event_driven`: no periodic tick. The model brings its own callbacks
  (`gnc_callbacks`). `BaseThrusterModel`, the burn scheduler, is event-driven.
"""
tick_schedule(model) = :periodic

# Scope: the spacecraft a guidance, navigation or control model acts on, and what
# happens to its command when one of them becomes inactive. Indices are run
# positions, from 1 in configuration order.

"""
    SpacecraftScope

The spacecraft a guidance, navigation or control model acts on, as its
`spacecraft_scope` declares: `AllSpacecraft()`, `BoundSpacecraft(indices)`,
`CoupledSpacecraft(primary, partners)` or a `SpacecraftSelection`.
"""
abstract type SpacecraftScope end

"""
    AllSpacecraft()

Scope of a model ticked for every spacecraft of the run, the default. Such a
model keeps one slot of run state per spacecraft.
"""
struct AllSpacecraft <: SpacecraftScope end

"""
    BoundSpacecraft(indices)

Scope of a model that acts only on the spacecraft at the run positions in the
tuple `indices`: the engine ticks it for those spacecraft only.
"""
struct BoundSpacecraft{N} <: SpacecraftScope
    indices::NTuple{N, Int}
end
BoundSpacecraft(indices::Tuple{Vararg{Integer}}) = BoundSpacecraft(map(Int, indices))

"""
    CoupledSpacecraft(primary, partners)

Scope of a model that couples several spacecraft: the engine ticks it for
`primary` only, and the tick also reads the spacecraft at the run positions in the
tuple `partners`, as the RPO controller ticks for the chaser and reads the target.
From SpaceAGORA 0.7.0 its output hooks and levels are called for `primary` only:
the partners are read, not acted on. A model that also acts on its partners
declares a scope that names them, such as `BoundSpacecraft` or `AllSpacecraft()`,
and is ticked for them too.
"""
struct CoupledSpacecraft{N} <: SpacecraftScope
    primary::Int
    partners::NTuple{N, Int}
end
CoupledSpacecraft(primary::Integer, partners::Tuple{Vararg{Integer}}) =
    CoupledSpacecraft(Int(primary), map(Int, partners))

"""
    spacecraft_scope(model) -> SpacecraftScope

The spacecraft guidance, navigation or control `model` acts on: every one, the
`BoundSpacecraft` indices, the `CoupledSpacecraft` primary, or the spacecraft a
`SpacecraftSelection` selects. Before the run the engine refuses a scope that names
a run position the run does not have. When the run starts, after the run-start
initializers, it takes the scope it reads then as the model's for the whole run; a
`SpacecraftSelection` is resolved then, its `where` called once per spacecraft.

During the run the engine calls the model only for the spacecraft its scope names:
its ticks and `release_command!`, its output hooks (`control_wrench`,
`control_propellant_rate`, `control_wheel_torque`), its levels (`actuator_levels`,
`control_thruster_levels`), `register_tstops!` and `atmosphere_transition!`.
SpaceAGORA 0.6 and earlier called the output hooks, the levels and
`atmosphere_transition!` for every spacecraft. Direct callers, such as the
conformance kit, may still call them for other spacecraft, so they keep returning
zeros or `nothing` there. `touchdown_spec` is still asked for every spacecraft
before the run, and so is `control_thruster_levels` when the result columns are
sized.

The scope is fixed once the run-start initializers ran. Each tick reads it again
(an event-driven model's scope is read once per accepted step), and a scope that
differs from the one the run took stops the run. So return the same scope at
every call: for a `SpacecraftSelection`, keep the selection in the model and return
that object, since its `where` function is compared by identity (`===`), and a
closure built anew at each call that captures new values is a different one. A dynamic effector's scope, read
once when the run starts, decides which spacecraft it acts on (see `scoped`).
Default: `AllSpacecraft()`.
"""
spacecraft_scope(model) = AllSpacecraft()

"""
    skip_inactive(model) -> Bool

Whether the engine stops ticking `model` for a spacecraft that is no longer active
(after an impact, a touchdown or another event that deactivates it). For a model
that returns `true`, while spacecraft `spacecraft_index` is inactive, each tick
time calls `release_command!(model, t, spacecraft_index)` instead of the tick. For
a `CoupledSpacecraft` model the primary decides; a model that must react to an
inactive partner checks `p.is_active` itself. Default: `true` (since SpaceAGORA
0.4.0; `false` before). Return `false` for a model that must keep ticking, as for
an active spacecraft.
"""
skip_inactive(model) = true

# Whether guidance `model` predicts its spacecraft's motion with the run's own
# gravity models, applied as to every spacecraft (EDG's prediction does). Private:
# the engine refuses to scope those gravity models in a run with such guidance
# (engine/effector_scopes.jl). A guidance type that predicts so adds the `true`
# method where it is defined, so the engine imports no module that loads after it.
_prediction_uses_run_gravity(model) = false

"""
    held_command(model, spacecraft_index) -> Any

The command control `model` holds for spacecraft `spacecraft_index` between its
ticks, which its output hooks read, as an immutable value or named tuple; the
model's own state is not shared with it. In `:mounted` actuation the levels the
mounted actuators apply, derived from it, are `actuator_levels(model,
spacecraft_index)`. Default: `nothing`, no held command, as for a spacecraft the
model does not act on.
"""
held_command(model, spacecraft_index) = nothing

"""
    held_command_time(model, spacecraft_index) -> Float64

The time (s into the run) of the tick that produced `held_command(model,
spacecraft_index)`. Default: `NaN`, before the first tick or for a model that does
not record it.
"""
held_command_time(model, spacecraft_index) = NaN

"""
    release_command!(model, t, spacecraft_index) -> nothing

Called at tick time `t` (s into the run) in place of the tick of a model whose
`skip_inactive` is `true`, while spacecraft `spacecraft_index` is inactive: set
the held command to the model's off command, such as zero thrust. Default:
nothing.
"""
release_command!(model, t, spacecraft_index) = nothing

"""
    gnc_callbacks(model, num_sats, args, period_s) -> Tuple

The callbacks an `:event_driven` control `model` adds to a run of `num_sats`
spacecraft with configuration `args`. `period_s` is the model's configured control
rate, for a model that also ticks periodically. The engine reads it only for
`:event_driven` models. Default: none.
"""
gnc_callbacks(model, num_sats, args, period_s) = ()

"""
    register_tstops!(integrator, model, spacecraft_index) -> nothing

Called after each periodic tick of control `model`, once for each spacecraft its
`spacecraft_scope` names, in ascending order (SpaceAGORA 0.5.1 and earlier called
it for every spacecraft index). Use it to add the solver stops that the new command
needs; the burn scheduler adds the ends of its burn window. The ticks of several
spacecraft may run on different threads, but these calls do not: they run one
spacecraft at a time, after all of that tick's spacecraft have ticked, because the
stop list belongs to the integrator they share. Default: none.
"""
register_tstops!(integrator, model, spacecraft_index) = nothing

# Links: the shared objects through which guidance and control models are paired,
# such as the RPO plan buffer or the energy-depletion state. A model declares the
# links it writes (`Produces`) and reads (`Consumes`). The engine groups the run's
# models by link identity, which survives the run's deep copy, and checks and
# initializes each link once.

"""
    Produces(link)

Wraps a shared object that a model writes for others to read, in `gnc_links`.
"""
struct Produces{L}
    link::L
end

"""
    Consumes(link)

Wraps a shared object that a model reads, in `gnc_links`. Before a periodic control
model's tick, the engine calls `before_reference_control!` for each model that
produces a link it consumes, unless the tick releases or skips every spacecraft.
"""
struct Consumes{L}
    link::L
end

"""
    gnc_links(model) -> Tuple

The links a guidance, navigation or control `model` takes part in, as `Produces`
and `Consumes` wrappers around the shared objects. Default: none.
"""
gnc_links(model) = ()

"""
    link_cardinality(link) -> NamedTuple

How many models may produce and consume `link`, as ranges:
`(producers = 0:1, consumers = 0:1)`. Before the run the engine refuses a link
whose counts fall outside them. Default: any number of each.
"""
link_cardinality(link) = (producers = 0:typemax(Int), consumers = 0:typemax(Int))

"""
    check_link(link, producers, consumers, args) -> nothing

Before the run, once per link: check that the models sharing `link` agree, such
as matching settings. Throw an `ArgumentError` to refuse the run. `producers` and
`consumers` are the configured models in configuration order. Default: nothing.
"""
check_link(link, producers, consumers, args) = nothing

"""
    link_initialize!(link, producers, consumers, u, p, t) -> nothing

At the start of the run, once per link and before any model's own initializer:
reset the run state the link owns. Default: nothing.
"""
link_initialize!(link, producers, consumers, u, p, t) = nothing

"""
    gnc_clock_started!(model, period_s, origin_s) -> nothing

Called with the period of control `model`'s periodic callback and the time its
clock starts from (s into the run): when the run starts, after the links are
initialized and before the models' own initializers, and again when the callback
starts, at every solve's start, so a checkpointed run calls it again at each
segment. The clock starts from the solve's start, or from the run's start (0) in a
checkpointed or resumed run, whose ticks keep one grid across segments. A model
that keeps a clock of its own ticks records it here. Default: nothing.
"""
gnc_clock_started!(model, period_s, origin_s) = nothing

"""
    gnc_period_s(model) -> Union{Nothing, Float64}

The tick period (s) that control `model`'s law assumes, such as the step its
controller was discretized with. Before the run the engine compares it with the
rate of the model's periodic callback and refuses the run when they differ.
Default: `nothing`, no declared period, as for a
law that measures the elapsed time itself.
"""
gnc_period_s(model) = nothing

end # module GNCRuntime
