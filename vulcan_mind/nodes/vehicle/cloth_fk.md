---
id: vehicle.cloth_fk
label: cloth_fk
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: cloth_fk
  lines:
  - 168
  - 209
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Robotics namespace supplying arm geometry, joint records, and kinematic
    conventions.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: pose
  type: ClothArmPose
  units: m,rad
  description: End-effector and link poses computed from the base pose and joint coordinates.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- robotics
- kinematics
charts:
- vehicle
origin: agent
---

# cloth_fk

## Purpose
`cloth_fk` computes the forward kinematics of the cloth-arm robotics model. It takes a base pose and joint vector, traverses the configured links and joints, and returns the pose data used by end-effector queries, inverse kinematics, and robot-arm dynamics.

## Theory & Math
For link transforms `T_i`, the end-effector transform is the ordered product `T_ee = T_base ∏ᵢ T_i(q_i)`. Joint rotations are applied about each configured joint axis, while link translations use the configured lengths and offsets. The output pose contains positions and orientations for the links and end effector.

## Model & Assumptions
The joint vector length must match the model’s joints, axes must be normalized or interpreted consistently by the rotation helper, and link geometry is treated as rigid. The base pose establishes the frame origin and orientation. Collision, actuator limits, and flexible deformation are outside this kinematic calculation.

## Design & Implementation
`robotics.jl` validates the joint vector, initializes the base transform, iterates through `ClothArmLink` and `ClothArmJoint` records, and accumulates each transform into a `ClothArmPose`. `cloth_fk_state` wraps this path for a `ClothArmState`; `cloth_ik` calls it repeatedly while evaluating trial configurations and the end-effector Jacobian.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Robotics namespace supplying arm geometry, joint records, and kinematic conventions. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `pose` | ClothArmPose | m,rad | — | End-effector and link poses computed from the base pose and joint coordinates. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_initial_state|cloth_robot_arm_initial_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:181-181`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:233-233`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_rest_quaternions|cloth_robot_arm_rest_quaternions]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:167-167`
- [[dynamics.cloth_robot_arm_dynamics_initialize_coupled_cloth_robot_arm_state_bang|initialize_coupled_cloth_robot_arm_state!]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:209-209`
- [[gnc.planner_core__robot_arm_plan_from_q_reference|_robot_arm_plan_from_q_reference]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:153-153`
- [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_base_wrench_ratios|_robot_arm_hypr_cloth_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:257-257`
- [[gnc.swarm_and_retiming__robot_arm_hypr_link_com_history|_robot_arm_hypr_link_com_history]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:156-156`
- [[gncz.clearance_robot_arm_clearance_stats_from_samples|robot_arm_clearance_stats_from_samples]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/clearance.jl:21-21`
- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:124-124`
- [[grp.src_gnc_robotics|gnc/robotics/]] · `members_out` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:153-153`
- [[vehicle.robotics__ee_position_jacobian|_ee_position_jacobian]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:251-251`
- [[vehicle.robotics_cloth_fk_state|cloth_fk_state]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:220-220`
- [[vehicle.robotics_cloth_ik|cloth_ik]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:274-274`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__quat_mul|_quat_mul]] · `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- `callees` → [[dynamics.cloth_multibody__rot|_rot]] · `callers` · call · `src/vehicle/robotics/robotics.jl:170-170`
- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_mul|_quat_mul]] · `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__rot|_rot]] · `callers` · call · `src/vehicle/robotics/robotics.jl:170-170`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/robotics/robotics.jl:189-189`
- `callees` → [[vehicle.robotics__quat_from_axis_angle|_quat_from_axis_angle]] · `callers` · call · `src/vehicle/robotics/robotics.jl:185-185`
- `callees` → [[vehicle.robotics__quat_mul|_quat_mul]] · `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- `callees` → [[vehicle.robotics__rot|_rot]] · `callers` · call · `src/vehicle/robotics/robotics.jl:170-170`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/vehicle/robotics/robotics.jl:186-186`
- `callees` → [[vehicle.robotics__validate_joint_vector|_validate_joint_vector]] · `callers` · call · `src/vehicle/robotics/robotics.jl:169-169`
- `callees` → [[vehicle.robotics_clotharmpose|ClothArmPose]] · `callers` · call · `src/vehicle/robotics/robotics.jl:198-198`
<!-- vulcan:connections:end -->

## Limitations
Forward kinematics reports the configured geometric model and does not certify physical reachability under joint stops or collision constraints. Large joint angles can expose convention errors if callers mix axis frames. Numerical orientation drift depends on the rotation representation and input values; downstream control should validate the returned pose before commanding an actuator.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl:168-209`.
