"""
    ExternalPropagation

Generic hook surface for spacecraft that an external integrator owns for a whole run (for example a
rigid-body or contact scene stepped at a fixed rate by an optional package). The engine keeps each owned
spacecraft in the state vector as a *shadow entry*: between syncs its right-hand side is the owner's chief
acceleration (translation) and plain attitude kinematics, and after every accepted solver step the engine
calls [`external_sync!`](@ref) and overwrites the shadow entries from the owner's absolute state.

This implementation module is engine-agnostic. Optional packages use its hook identities through the
versioned `SpaceAGORA.ExternalPropagationServices` facade, subtype [`AbstractExternalPropagator`](@ref),
and attach owners with that facade's `with_external_propagators`. Defaults throw `Not implemented` until a
method is added, except for [`external_preflight`](@ref), whose default accepts.

Contract, per run:

1. `external_preflight(ep, args)` may throw `ArgumentError` for owner-specific problems (called once at setup,
   after the engine's own refusals).
2. `external_prepare(ep, args, u0)` returns the run's mutable *runtime*. It must not alias `ep` or any state
   another run can reach (a campaign may share `ep` across threads), and it is built from the initial state
   `u0`, which is the single source of truth for the owned spacecraft's initial conditions.
3. `external_sync!(rt, t)` advances the owner to time `t` or the last owned step not past it, and returns the
   number of owner steps taken. The engine calls it after every accepted step, before any other discrete
   callback, so guidance, navigation and control callbacks see state current to `t` at an owner-step boundary.
4. `external_state(rt, t, slot)` returns the absolute state of the owned spacecraft at position `slot` of the
   owner's `external_spacecraft` list, at engine time `t`, in SpaceAGORA's conventions: `pos`, `vel` (inertial,
   m and m/s), `q` (scalar last, `v_body = rot(q) * v_inertial`) and `ω` (body frame). The owner's time is at
   most one owner step behind `t`; the owner advances its state over that remainder with the same model the
   shadow right-hand side uses (chief acceleration, constant body rate), so a written shadow entry is consistent
   in time with the state the solver integrates.
5. `external_acceleration(rt, t, slot)` is the translational acceleration the shadow entry follows at engine time
   `t` until the next sync (inertial, m/s^2); an owner uses its chief's acceleration, evaluated at `t`.

Contract 2.0.0 puts the time before the owner-local position `slot`, as every SpaceAGORA hook puts time before
an index; contract 1.0.0 had `external_state(rt, k, t)` and `external_acceleration(rt, k, t)`.
"""
module ExternalPropagation

using StaticArrays

export AbstractExternalPropagator
export external_spacecraft, external_step, external_preflight, external_prepare
export external_sync!, external_time, external_state, external_acceleration
export external_capabilities, external_stage_safe, external_bind_force_service!, external_control_changed!
export external_checkpoint_key, external_checkpoint_state, external_restore_state!
export external_stage_key


"""
    AbstractExternalPropagator

Supertype for a template that owns selected spacecraft throughout a simulation.
Use the versioned `SpaceAGORA.ExternalPropagationServices` facade to implement
the eight `external_*` hooks. The engine keeps owned spacecraft as shadow entries
for ordinary output and guidance, with owner-supplied translational acceleration
and attitude kinematics between synchronization points.

Templates may be shared by concurrent runs. `external_prepare` must construct a
separate mutable runtime from the supplied initial state; neither preparation
nor stepping may mutate the template or another run's runtime. After every
accepted solver step, `external_sync!` advances whole owner steps without passing
the solver time, before other discrete callbacks run. `external_state` supplies
the state at solver time, including any remainder shorter than one owner step,
using a model consistent with `external_acceleration` and attitude kinematics.

Ownership, solver, callback-rate and event limits are checked by the facade's
`validate_external_propagators` and by run setup. No optional engine is loaded
merely by loading SpaceAGORA or this interface.
"""
abstract type AbstractExternalPropagator end

_unimplemented(f, T) = throw(ErrorException("Not implemented: $(f) for $(T); implement it for the AbstractExternalPropagator subtype."))

"""Opt in to experimental lifecycle services; legacy owners retain all existing refusals.
`configured_forces` requires the translation-only force-service contract; `checkpoint` requires
versioned owner snapshots. Neither flag asserts numerical validity of an adapter.
"""
external_capabilities(::AbstractExternalPropagator) = (configured_forces=false, checkpoint=false)

"""Explicit declaration that a model's force hooks are pure trial-state evaluations, without
controller ticks, history advancement, mass loss or torque. Unknown models are refused by default.
The initial force service also admits the built-in central and J2 gravity models.
"""
external_stage_safe(model) = false

"""Stable configuration key for a declared stage-safe custom force/control model. Include all
parameters affecting force or control, excluding held commands and other checkpointed runtime state.
Required for configured-force owners even before checkpointing is enabled.
"""
external_stage_key(model) = _unimplemented("external_stage_key", typeof(model))

"""Attach this run's opaque callable force service to an opted-in runtime. Call `service(pos,vel,t)`
for inertial acceleration in SI at physical elapsed time; it includes configured held control.
Do not retain the service across runs or serialize it. Trial evaluation must not advance the owner.
"""
external_bind_force_service!(rt, service) = _unimplemented("external_bind_force_service!", typeof(rt))

"""Invalidate owner force caches after a physical control tick. The owner has already reached `t`
using the old command; the next interval must use the newly held command. No stepping is allowed here.
"""
external_control_changed!(rt, t::Float64) = _unimplemented("external_control_changed!", typeof(rt))

"""Stable, serializable adapter/schema/configuration identity for an opted-in owner snapshot.
Include adapter version and every owner setting that affects its trajectory; do not use runtime addresses.
"""
external_checkpoint_key(ep::AbstractExternalPropagator) = _unimplemented("external_checkpoint_key", typeof(ep))

"""Return an independent serializable snapshot of the owner's mathematical state at physical `t`.
Include native coordinates, clock/grid origin and required numerical state, but not the force service.
"""
external_checkpoint_state(rt, t::Float64) = _unimplemented("external_checkpoint_state", typeof(rt))

"""Restore a validated owner snapshot at physical `t`, after the engine restores model state and
before any resumed propagation or output. Preserve this run's force service and do not tick controls.
"""
external_restore_state!(rt, state, t::Float64) = _unimplemented("external_restore_state!", typeof(rt))


"""`external_spacecraft(ep) -> Vector{Int}`: run indices (1-based, in `dynamics_model.spacecraft`) the owner propagates."""
external_spacecraft(ep::AbstractExternalPropagator) = _unimplemented("external_spacecraft", typeof(ep))

"""`external_step(ep) -> Float64`: the owner's fixed step [s]. GNC rates must be integer multiples of it."""
external_step(ep::AbstractExternalPropagator) = _unimplemented("external_step", typeof(ep))

"""`external_preflight(ep, args)`: owner-specific setup checks (throw `ArgumentError`); the default accepts."""
external_preflight(::AbstractExternalPropagator, args) = nothing

"""
    external_prepare(ep, args, u0) -> runtime

Build the run's own mutable runtime from its initial state. Read spacecraft
initial conditions through `ExternalPropagationServices.external_initial_state`;
the engine's storage layout is not part of the adapter contract. The result must
not alias mutable state in the owner template or another run.
"""
external_prepare(ep::AbstractExternalPropagator, args, u0) = _unimplemented("external_prepare", typeof(ep))

"""`external_sync!(rt, t) -> Int`: advance by whole owner steps to the last step at or before elapsed time `t` [s]; return the number of steps taken."""
external_sync!(rt, t::Float64) = _unimplemented("external_sync!", typeof(rt))

"""`external_time(rt) -> Float64`: the owner's current time [s], an integer number of steps from its epoch. It lets tests and tools observe the owner's clock. The engine calls it only for an owner that opts into the experimental lifecycle extension (`external_capabilities`), to check that the owner reached a control tick or checkpoint time."""
external_time(rt) = _unimplemented("external_time", typeof(rt))

"""
    external_state(rt, t, slot) -> (pos, vel, q, ω)

Return the absolute state, at elapsed engine time `t` [s], of the owned spacecraft
at position `slot` (1-based) of the owner's `external_spacecraft` list; `slot` is
local to the owner, not a run index. The fields are inertial position `pos` [m],
inertial velocity `vel` [m/s], attitude quaternion `q` and body angular velocity
`ω` [rad/s]. Account for any lag after the last whole owner step with the model
used by the shadow acceleration and attitude kinematics.

`q = (x, y, z, w)` is scalar last and gives `v_body = rot(q) * v_inertial`, with
`rot(q) = (w^2 - e⋅e) I + 2 e eᵀ - 2 w [e×]`, `e = (x, y, z)` and `[e×]` the
cross-product matrix of `e`. This is the engine's convention for the state `q`,
`StateSample.q_ib` and the saved `sc<i>_q_*` columns.

Contract 2.0.0 order: time before the position. Contract 1.0.0 was
`external_state(rt, k, t)`.
"""
external_state(rt, t::Float64, slot::Int) = _unimplemented("external_state", typeof(rt))

"""
    external_acceleration(rt, t, slot) -> SVector{3,Float64}

Inertial translational acceleration [m/s^2] that the shadow entry of the owned
spacecraft at position `slot` follows at elapsed engine time `t` [s]. Contract
2.0.0 order: time before the position; contract 1.0.0 was
`external_acceleration(rt, k, t)`.
"""
external_acceleration(rt, t::Float64, slot::Int) = _unimplemented("external_acceleration", typeof(rt))

end # module ExternalPropagation
