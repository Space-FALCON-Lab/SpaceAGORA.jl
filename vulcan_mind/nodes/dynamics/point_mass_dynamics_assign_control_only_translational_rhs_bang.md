---
id: dynamics.point_mass_dynamics_assign_control_only_translational_rhs_bang
label: assign_control_only_translational_rhs!
kind: function
source:
  file: src/dynamics/translational/point_mass_dynamics.jl
  symbol: assign_control_only_translational_rhs!
  lines:
  - 45
  - 45
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
- id: mass_rate
  type: Real
  units: n/a
  required: true
  description: Positional argument `mass_rate`.
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
  description: Return value of `assign_control_only_translational_rhs!`; mutates `du_view`
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

# assign_control_only_translational_rhs!

## Purpose

`assign_control_only_translational_rhs!` writes a translational derivative in which position is frozen but velocity and mass respond to control forces. It supports the control-only channel of a split integrator, where the kinematic position update is applied elsewhere.

## Design & Implementation

The position slot is filled from `zero_position_derivative()` rather than from `sc_view.vel`, so `du_view.pos` receives an explicit zero vector. Velocity comes from `acceleration_from_force(net_force, sc_view.mass)` and mass from `mass_derivative(mass_rate)`, matching the full variant. The function is `@inline`, mutates `du_view` in place and returns `nothing`; the zero-position helper keeps the assignment allocation-free instead of building a fresh zero vector at each call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du_view` | Any | n/a | yes | Positional argument `du_view`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `net_force` | AbstractVector{<:Real} | n/a | yes | Positional argument `net_force`. |
| in | `mass_rate` | Real | n/a | yes | Positional argument `mass_rate`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `assign_control_only_translational_rhs!`; mutates `du_view` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- [[simulation.dynamics_rhs_spacecraft_dynamics_fast_control_bang|spacecraft_dynamics_fast_control!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2223-2223`

**Downstream**

- `callees` → [[dynamics.point_mass_dynamics_mass_derivative|mass_derivative]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:53-53`
- `callees` → [[dynamics.position_kinematics_zero_position_derivative|zero_position_derivative]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:51-51`
- `callees` → [[dynx.translational_point_mass_dynamics_acceleration_from_force|acceleration_from_force]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:52-52`
<!-- vulcan:connections:end -->

## Limitations

Zeroing the position derivative is only physically meaningful inside a splitting scheme that adds the kinematic term on another pass; invoked as the sole right-hand side it integrates a vehicle that accelerates but never moves. As with the other variants, `du_view` and `sc_view` are untyped, and a non-finite or near-zero `sc_view.mass` yields a silently zeroed acceleration through the guard inside `acceleration_from_force`.

## Provenance
Mapped from `src/dynamics/translational/point_mass_dynamics.jl` line 45.
