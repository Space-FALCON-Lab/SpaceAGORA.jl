"""
    ExternalPropagation

Generic hook surface for spacecraft that an external integrator owns for a whole run (for example a
rigid-body or contact scene stepped at a fixed rate by an optional package). The engine keeps each owned
spacecraft in the state vector as a *shadow entry*: between syncs its right-hand side is the owner's chief
acceleration (translation) and plain attitude kinematics, and after every accepted solver step the engine
calls [`external_sync!`](@ref) and overwrites the shadow entries from the owner's absolute state.

This module is MuJoCo-agnostic. It is an internal `SimulationModel` name space (not exported from the root
package); an implementation subtypes [`AbstractExternalPropagator`](@ref) and is passed in
`SimulationConfiguration.external_propagators`. Nothing here defines behavior for any concrete owner, so every
generic below throws `Not implemented` until a method is added.

Contract, per run:

1. `external_preflight(ep, args)` may throw `ArgumentError` for owner-specific problems (called once at setup,
   after the engine's own refusals).
2. `external_prepare(ep, args, u0)` returns the run's mutable *runtime*. It must not alias `ep` or any state
   another run can reach (a campaign may share `ep` across threads), and it is built from the initial state
   `u0`, which is the single source of truth for the owned spacecraft's initial conditions.
3. `external_sync!(rt, t)` advances the owner to time `t` or the last owned step not past it, and returns the
   number of owner steps taken. The engine calls it after every accepted step, before any other discrete
   callback, so guidance, navigation and control callbacks see state current to `t` at an owner-step boundary.
4. `external_state(rt, k, t)` returns the absolute state of the `k`-th owned spacecraft at engine time `t` in
   SpaceAGORA's conventions: `pos`, `vel` (inertial, m and m/s), `q` (scalar-last, body to inertial) and `ω`
   (body frame). The owner's time is at most one owner step behind `t`; the owner advances its state over that
   remainder with the same model the shadow right-hand side uses (chief acceleration, constant body rate), so
   a written shadow entry is consistent in time with the state the solver integrates.
5. `external_acceleration(rt, k, t)` is the translational acceleration the shadow entry follows at engine time
   `t` until the next sync (inertial, m/s^2); an owner uses its chief's acceleration, evaluated at `t`.
"""
module ExternalPropagation

using StaticArrays

export AbstractExternalPropagator
export external_spacecraft, external_step, external_preflight, external_prepare
export external_sync!, external_time, external_state, external_acceleration

"""Supertype of externally propagated spacecraft owners; see [`ExternalPropagation`](@ref)."""
abstract type AbstractExternalPropagator end

_unimplemented(f, T) = throw(ErrorException("Not implemented: $(f) for $(T); implement it for the AbstractExternalPropagator subtype."))

"""`external_spacecraft(ep) -> Vector{Int}`: run indices (1-based, in `dynamics_model.spacecraft`) the owner propagates."""
external_spacecraft(ep::AbstractExternalPropagator) = _unimplemented("external_spacecraft", typeof(ep))

"""`external_step(ep) -> Float64`: the owner's fixed step [s]. GNC rates must be integer multiples of it."""
external_step(ep::AbstractExternalPropagator) = _unimplemented("external_step", typeof(ep))

"""`external_preflight(ep, args)`: owner-specific setup checks (throw `ArgumentError`); the default accepts."""
external_preflight(::AbstractExternalPropagator, args) = nothing

"""`external_prepare(ep, args, u0) -> runtime`: build the run's own mutable runtime from the initial state."""
external_prepare(ep::AbstractExternalPropagator, args, u0) = _unimplemented("external_prepare", typeof(ep))

"""`external_sync!(rt, t) -> Int`: advance the owner to time `t` (elapsed seconds); returns the owner steps taken."""
external_sync!(rt, t::Float64) = _unimplemented("external_sync!", typeof(rt))

"""`external_time(rt) -> Float64`: the owner's current time [s], an integer number of steps from its epoch."""
external_time(rt) = _unimplemented("external_time", typeof(rt))

"""`external_state(rt, k, t) -> (pos, vel, q, ω)`: absolute state of the `k`-th owned spacecraft at engine time `t`, SpaceAGORA conventions."""
external_state(rt, k::Int, t::Float64) = _unimplemented("external_state", typeof(rt))

"""`external_acceleration(rt, k, t) -> SVector{3,Float64}`: shadow-entry translational acceleration at engine time `t`."""
external_acceleration(rt, k::Int, t::Float64) = _unimplemented("external_acceleration", typeof(rt))

end # module ExternalPropagation
