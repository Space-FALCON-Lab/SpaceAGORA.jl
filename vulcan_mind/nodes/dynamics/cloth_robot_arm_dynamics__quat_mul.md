---
id: dynamics.cloth_robot_arm_dynamics__quat_mul
label: _quat_mul
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _quat_mul
  lines:
  - 61
  - 61
inputs:
- id: a
  type: Any
  units: n/a
  required: true
  description: Positional argument `a`.
- id: b
  type: Any
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
- dynamics
charts:
- dynamics
origin: agent
---

# _quat_mul

## Purpose
Composes two rotations by Hamilton product in the module's scalar-last convention, normalising inputs and output.

## Theory & Math
For $a=(\mathbf{a}_v, a_w)$ and $b=(\mathbf{b}_v, b_w)$ with vector parts $\mathbf{a}_v,\mathbf{b}_v\in\mathbb{R}^3$ and scalars $a_w,b_w$, the product is $ab = (a_w\mathbf{b}_v + b_w\mathbf{a}_v + \mathbf{a}_v\times\mathbf{b}_v,\; a_w b_w - \mathbf{a}_v\cdot\mathbf{b}_v)$, which the code then divides by its norm.

## Design & Implementation
Both operands pass through `_unit_quat`, are destructured as `(x, y, z, w)`, and the product components are `aw*bx + bw*ax + ay*bz - az*by`, `aw*by + bw*ay + az*bx - ax*bz`, `aw*bz + bw*az + ax*by - ay*bx`, and `aw*bw - ax*bx - ay*by - az*bz`. The result is normalised again before return. This is the same formula as `_quat_raw_mul` but with three normalisations.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Any | n/a | yes | Positional argument `a`. |
| in | `b` | Any | n/a | yes | Positional argument `b`. |
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
- [[grp.src_vehicle_robotics|vehicle/robotics/]] · `members_out` → `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:62-62`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:62-62`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:62-62`
<!-- vulcan:connections:end -->

## Limitations
Three normalisations per multiply add roughly 3 square roots and divisions, noticeable inside the per-link RHS loop. Normalising the output makes the function unsuitable for quaternion derivatives (hence the separate `_quat_raw_mul`). Degenerate inputs are replaced by identity without error.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 61.
