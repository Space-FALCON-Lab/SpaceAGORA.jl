---
id: dynamics.cloth_robot_arm_dynamics__rot
label: _rot
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _rot
  lines:
  - 75
  - 75
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
- dynamics
charts:
- dynamics
origin: agent
---

# _rot

## Purpose
Converts a scalar-last quaternion into the 3 x 3 body-to-world rotation matrix used throughout the coupled arm dynamics.

## Theory & Math
$R(q) = \begin{pmatrix} 1-2(y^2+z^2) & 2(xy - zw) & 2(xz + yw)\\ 2(xy + zw) & 1-2(x^2+z^2) & 2(yz - xw)\\ 2(xz - yw) & 2(yz + xw) & 1-2(x^2+y^2)\end{pmatrix}$ for unit $q=(x,y,z,w)$.

## Design & Implementation
Normalises with `_unit_quat`, destructures `(x, y, z, w)`, and builds an `SMatrix{3,3,Float64}` in column-major order. The first column is `(1-2y²-2z², 2xy+2zw, 2xz-2yw)`, the second `(2xy-2zw, 1-2x²-2z², 2yz+2xw)`, and the third `(2xz+2yw, 2yz-2xw, 1-2x²-2y²)`. Callers use `R * v_body` to map into the world frame and `R' * v_world` for the inverse.

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
- [[grp.src_vehicle_robotics|vehicle/robotics/]] · `members_out` → `callers` · call · `src/vehicle/robotics/robotics.jl:170-170`
- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:170-170`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:76-76`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:76-76`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:76-76`
<!-- vulcan:connections:end -->

## Limitations
The normalisation cost is paid on every call, and `_rot` is called several times per link per RHS evaluation. A degenerate quaternion yields the identity matrix silently. The matrix is only orthonormal to the precision of the normalised input.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 75.
