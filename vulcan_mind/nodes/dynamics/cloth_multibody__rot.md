---
id: dynamics.cloth_multibody__rot
label: _rot
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _rot
  lines:
  - 154
  - 154
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
Converts a quaternion to the body-to-world rotation matrix.

## Design & Implementation
Normalises, destructures into `x, y, z, w`, and builds the standard `SMatrix{3,3}` in column-major order. `@inline`. Called for every body and attachment point in each derivative evaluation.

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

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:155-155`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:155-155`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:155-155`
<!-- vulcan:connections:end -->

## Limitations
Recomputed rather than cached per body per evaluation; `_parent_kinematics` is called twice per joint so each body's rotation is rebuilt several times per step.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 154.
