---
id: vehicle.robotics__quat_mul
label: _quat_mul
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: _quat_mul
  lines:
  - 82
  - 82
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
  description: Return value of `_quat_mul`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# _quat_mul

## Purpose
Composes two scalar-last quaternions in the project's convention, used to accumulate link orientation down the chain.

## Theory & Math
For $a = (a_v, a_w)$ and $b = (b_v, b_w)$:

$$
(ab)_v = a_w b_v + b_w a_v + a_v \times b_v,\qquad (ab)_w = a_w b_w - a_v \cdot b_v
$$

## Design & Implementation
Destructures both inputs into `x, y, z, w` and returns the Hamilton product with the vector part first and scalar last, written out as four explicit expressions so there is no temporary array. `@inline` on static vectors.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | SVector{4, Float64} | n/a | yes | Positional argument `a`. |
| in | `b` | SVector{4, Float64} | n/a | yes | Positional argument `b`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{4, | n/a | — | Return value of `_quat_mul`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody__rest_child_parent_quat|_rest_child_parent_quat]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:547-547`
- [[dynamics.cloth_robot_arm_dynamics__rest_child_parent_quat|_rest_child_parent_quat]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:161-161`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:392-392`
- [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_out` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:277-277`
- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The product is not re-normalised, so repeated composition drifts from unit length; `cloth_fk` wraps the result in `_unit_quat` at each link for that reason.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 82.
