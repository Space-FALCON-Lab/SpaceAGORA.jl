---
id: dynamics.cloth_robot_arm_dynamics__quat_raw_mul
label: _quat_raw_mul
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _quat_raw_mul
  lines:
  - 148
  - 148
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
Hamilton product of two `SVector{4,Float64}` quaternions without normalisation, required for computing the quaternion time derivative.

## Theory & Math
$\dot q = \tfrac{1}{2}\, q \otimes (\boldsymbol\omega_b, 0)$, where $q$ is the body-to-world attitude quaternion (scalar last) and $\boldsymbol\omega_b$ is the angular velocity in the body frame.

## Design & Implementation
Identical component formulas to `_quat_mul` but typed strictly on `SVector{4,Float64}` and with no `_unit_quat` calls. `assign_coupled_cloth_robot_arm_rhs!` uses it as `0.5 .* _quat_raw_mul(q, (ω_x, ω_y, ω_z, 0))` to form `q̇` from the body angular velocity.

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
No check that inputs are unit, so the derivative inherits any drift in `q`; the integrator is expected to renormalise through `_unit_quat` on reads. Because the derivative uses the raw `q` from the state, a state quaternion with norm far from 1 gives a scaled derivative.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 148.
