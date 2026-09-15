---
id: dynamics.cloth_multibody__quat_raw_mul
label: _quat_raw_mul
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _quat_raw_mul
  lines:
  - 142
  - 142
inputs:
- id: a
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Positional argument `a`.
- id: b
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Positional argument `b`.
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
  type: SVector{4,
  units: n/a
  description: Return value of `_quat_raw_mul`.
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

# _quat_raw_mul

## Purpose
Composes two quaternions without normalising, required for computing the quaternion derivative where the pure-vector rate factor is not a unit quaternion.

## Theory & Math
$$
\dot{q} = \tfrac{1}{2}\, q \otimes (\omega_b, 0)
$$

## Design & Implementation
Destructures both `SVector{4}` inputs and returns the Hamilton product directly. `@inline`. Used in `compliant_multibody_dynamics` as `0.5 * q ⊗ (ω, 0)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | SVector{4, Float64} | n/a | yes | Positional argument `a`. |
| in | `b` | SVector{4, Float64} | n/a | yes | Positional argument `b`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{4, | n/a | — | Return value of `_quat_raw_mul`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:417-417`
- [[dynx.multibody_cloth_compliant_multibody_dynamics|compliant_multibody_dynamics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:618-618`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Typed strictly on `SVector{4,Float64}`, so callers must convert first.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 142.
