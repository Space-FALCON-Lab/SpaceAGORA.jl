---
id: dynamics.position_kinematics_zero_position_derivative
label: zero_position_derivative
kind: function
source:
  file: src/dynamics/translational/position_kinematics.jl
  symbol: zero_position_derivative
  lines:
  - 7
  - 7
inputs:
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
  type: SVector{3,
  units: n/a
  description: Return value of `zero_position_derivative`.
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

# zero_position_derivative

## Purpose

Returns the zero translational velocity vector used as the position derivative when a satellite's orbital state is frozen or when a propagator needs a neutral element for the position block of the state derivative.

## Design & Implementation

A one-line `@inline` function taking no arguments and returning `SVector{3, Float64}(0.0, 0.0, 0.0)`. Because the return type is a statically sized `SVector`, the constant is stack allocated and folds away at compile time, so calling it inside a derivative kernel costs nothing. It is the companion of `position_derivative`, which converts an arbitrary velocity vector into the same `SVector{3, Float64}` type; using one shared return type keeps the derivative assembly type stable.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `zero_position_derivative`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.point_mass_dynamics_assign_control_only_translational_rhs_bang|assign_control_only_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:51-51`
- [[dynamics.point_mass_dynamics_assign_force_only_translational_rhs_bang|assign_force_only_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:62-62`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/translational/position_kinematics.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Fixed at three components and `Float64`; it cannot serve dual-number or `BigFloat` derivative pipelines used for automatic differentiation or extended-precision propagation. It carries no notion of reference frame or units, so the caller must guarantee the zero it splices into the state derivative belongs in the same frame (metres per second, inertial) as the surrounding position rows.

## Provenance
Mapped from `src/dynamics/translational/position_kinematics.jl` line 7.
