---
id: dynamics.point_mass_dynamics_assign_slow_translational_rhs_bang
label: assign_slow_translational_rhs!
kind: function
source:
  file: src/dynamics/translational/point_mass_dynamics.jl
  symbol: assign_slow_translational_rhs!
  lines:
  - 34
  - 34
inputs:
- id: du_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `du_view`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: net_force
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `net_force`.
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: Nothing
  units: n/a
  description: Return value of `assign_slow_translational_rhs!`; mutates `du_view`
    in place. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# assign_slow_translational_rhs!

## Purpose

`assign_slow_translational_rhs!` fills the translational derivative for the slow half of a multi-rate integration split: position and velocity still evolve, but the mass derivative is pinned to zero because propellant flow is handled on the fast channel.

## Design & Implementation

It broadcasts `position_derivative(sc_view.vel)` into `du_view.pos` and `acceleration_from_force(net_force, sc_view.mass)` into `du_view.vel`, exactly as the full variant does, then assigns the literal `du_view.mass = 0.0` instead of calling `mass_derivative`. There is no `mass_rate` argument in the signature at all, which makes the omission explicit at every call site. The function is `@inline`, mutates `du_view`, and returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du_view` | Any | n/a | yes | Positional argument `du_view`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `net_force` | AbstractVector{<:Real} | n/a | yes | Positional argument `net_force`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `assign_slow_translational_rhs!`; mutates `du_view` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1378-1378`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1869-1869`

**Downstream**

- `callees` → [[dynx.translational_point_mass_dynamics_acceleration_from_force|acceleration_from_force]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:40-40`
- `callees` → [[dynx.translational_position_kinematics_position_derivative|position_derivative]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:39-39`
<!-- vulcan:connections:end -->

## Limitations

Holding `du_view.mass` at zero is only correct when some other right-hand side owns the mass channel for the same step; used on its own it produces a vehicle that burns propellant without losing mass. Because `sc_view.mass` is still read for the acceleration, the two channels must stay synchronised or the velocity derivative uses a stale mass. The `du_view` and `sc_view` arguments are untyped, so field-name mismatches only fail at runtime.

## Provenance
Mapped from `src/dynamics/translational/point_mass_dynamics.jl` line 34.
