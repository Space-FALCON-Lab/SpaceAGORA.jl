---
id: dynamics.cloth_robot_arm_dynamics__unit_quat
label: _unit_quat
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: _unit_quat
  lines:
  - 48
  - 48
inputs:
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
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
  description: Return value of `_unit_quat`.
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

# _unit_quat

## Purpose
Normalises any 4-element quaternion-like input to a unit `SVector{4,Float64}`, substituting the identity quaternion when the input is degenerate.

## Design & Implementation
Converts `q` to `SVector{4,Float64}`, computes `nq = norm(qv)`, and returns `qv / nq` when `nq` is finite and greater than `eps(Float64)`; otherwise returns `(0,0,0,1)`. Marked `@inline`; used by every quaternion routine in the module to guard against drift accumulated by the integrator.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{4, | n/a | — | Return value of `_unit_quat`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody__axis_angle_error|_axis_angle_error]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:166-166`
- [[dynamics.cloth_multibody__joint_rest|_joint_rest]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:511-511`
- [[dynamics.cloth_multibody__normalize_state_quaternions__normalize_state_quaternions_bang|_normalize_state_quaternions!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:473-473`
- [[dynamics.cloth_multibody__quat_conj|_quat_conj]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:123-123`
- [[dynamics.cloth_multibody__quat_mul|_quat_mul]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:129-129`
- [[dynamics.cloth_multibody__rot|_rot]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:155-155`
- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:289-289`
- [[dynamics.cloth_multibody_compliant_state_parts|compliant_state_parts]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:462-462`
- [[dynamics.cloth_multibody_compliant_state_vector|compliant_state_vector]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:450-450`
- [[dynamics.cloth_multibody_rectangular_prism_inertia|rectangular_prism_inertia]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:217-217`
- [[dynamics.cloth_robot_arm_dynamics__axis_angle_error|_axis_angle_error]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:136-136`
- [[dynamics.cloth_robot_arm_dynamics__coupled_body_state|_coupled_body_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:299-299`
- [[dynamics.cloth_robot_arm_dynamics__coupled_parent_kinematics|_coupled_parent_kinematics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:308-308`
- [[dynamics.cloth_robot_arm_dynamics__quat_conj|_quat_conj]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:56-56`
- [[dynamics.cloth_robot_arm_dynamics__quat_mul|_quat_mul]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:62-62`
- [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:76-76`
- [[dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang|initialize_coupled_cloth_robot_arm_state!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:205-205`
- [[dynx.multibody_cloth_assign_coupled_cloth_robot_arm_rhs_bang|assign_coupled_cloth_robot_arm_rhs!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:360-360`
- [[grp.src_vehicle_robotics|vehicle/robotics/]] · `members_out` → `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- [[vehicle.robotics__rot|_rot]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:106-106`
- [[vehicle.robotics_clotharmbasepose|ClothArmBasePose]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:26-26`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Silently replacing a zero or NaN quaternion with identity hides state corruption rather than raising an error. The threshold `eps(Float64)` (about 2.2e-16) accepts extremely small but nonzero vectors, whose normalisation amplifies noise. Sign is preserved, so `q` and `-q` remain distinct even though they represent the same rotation.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 48.
