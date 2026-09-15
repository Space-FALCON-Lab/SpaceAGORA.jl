---
id: vehicle.robotics__rot
label: _rot
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: _rot
  lines:
  - 105
  - 105
inputs:
- id: q_raw
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_raw`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `_rot`.
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

# _rot

## Purpose
Converts a quaternion to the three-by-three rotation matrix that maps parent-frame vectors into the world frame.

## Theory & Math
$$
R = \begin{bmatrix} 1-2(y^2+z^2) & 2(xy - zw) & 2(xz + yw) \\ 2(xy + zw) & 1-2(x^2+z^2) & 2(yz - xw) \\ 2(xz - yw) & 2(yz + xw) & 1-2(x^2+y^2) \end{bmatrix}
$$

## Design & Implementation
Normalises the input through `_unit_quat`, destructures it, and constructs an `SMatrix{3,3}` from the standard quadratic form in `x, y, z, w`. The literal is written column-major to match `SMatrix` construction order. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_raw` | Any | n/a | yes | Positional argument `q_raw`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_rot`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody__body_offset_velocity|_body_offset_velocity]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:179-179`
- [[dynamics.cloth_multibody__parent_kinematics|_parent_kinematics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:485-485`
- [[dynamics.cloth_robot_arm_dynamics__body_offset_velocity_world|_body_offset_velocity_world]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:131-131`
- [[dynamics.cloth_robot_arm_dynamics__coupled_parent_kinematics|_coupled_parent_kinematics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:312-312`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_end_effector|cloth_robot_arm_end_effector]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:433-433`
- [[dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang|initialize_coupled_cloth_robot_arm_state!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:210-210`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:361-361`
- [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_out` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:179-179`
- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:170-170`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:106-106`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:106-106`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:106-106`
<!-- vulcan:connections:end -->

## Limitations
Normalising on every call costs a square root even when the caller already holds a unit quaternion; FK calls this at least twice per link.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 105.
