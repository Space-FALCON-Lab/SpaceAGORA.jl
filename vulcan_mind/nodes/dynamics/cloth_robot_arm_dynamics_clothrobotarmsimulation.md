---
id: dynamics.cloth_robot_arm_dynamics_clothrobotarmsimulation
label: ClothRobotArmSimulation
kind: struct
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: ClothRobotArmSimulation
  lines:
  - 36
  - 36
inputs:
- id: multibody_model
  type: CompliantMultibodyModel
  units: n/a
  required: true
  description: Field `multibody_model`.
- id: trajectory
  type: CompliantMultibodyTrajectory
  units: n/a
  required: true
  description: Field `trajectory`.
- id: end_effector_positions
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `end_effector_positions`.
- id: reference_end_effector_positions
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `reference_end_effector_positions`.
- id: tracking_error_m
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `tracking_error_m`.
- id: joint_compliance_torques_body
  type: Array{Float64, 3}
  units: n/a
  required: true
  description: Field `joint_compliance_torques_body`.
- id: joint_actuator_torques_body
  type: Array{Float64, 3}
  units: n/a
  required: true
  description: Field `joint_actuator_torques_body`.
- id: joint_total_torques_body
  type: Array{Float64, 3}
  units: n/a
  required: true
  description: Field `joint_total_torques_body`.
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
  type: ClothRobotArmSimulation
  units: n/a
  description: Constructed `ClothRobotArmSimulation`.
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

# ClothRobotArmSimulation

## Purpose
Result container from `simulate_cloth_robot_arm_plan`, bundling the compliant model, the state trajectory, end-effector tracking data, and per-joint torque histories.

## Design & Implementation
Fields: `multibody_model::CompliantMultibodyModel`, `trajectory::CompliantMultibodyTrajectory` (times and packed state vectors), `end_effector_positions::Matrix{Float64}` and `reference_end_effector_positions::Matrix{Float64}` (both 3 x N_t, metres), `tracking_error_m::Vector{Float64}` (Euclidean norm per sample), and three `Array{Float64,3}` of shape 3 x n_joints x N_t holding compliance, actuator, and total joint torques expressed in the child body frame (N m).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `multibody_model` | CompliantMultibodyModel | n/a | yes | Field `multibody_model`. |
| in | `trajectory` | CompliantMultibodyTrajectory | n/a | yes | Field `trajectory`. |
| in | `end_effector_positions` | Matrix{Float64} | n/a | yes | Field `end_effector_positions`. |
| in | `reference_end_effector_positions` | Matrix{Float64} | n/a | yes | Field `reference_end_effector_positions`. |
| in | `tracking_error_m` | Vector{Float64} | n/a | yes | Field `tracking_error_m`. |
| in | `joint_compliance_torques_body` | Array{Float64, 3} | n/a | yes | Field `joint_compliance_torques_body`. |
| in | `joint_actuator_torques_body` | Array{Float64, 3} | n/a | yes | Field `joint_actuator_torques_body`. |
| in | `joint_total_torques_body` | Array{Float64, 3} | n/a | yes | Field `joint_total_torques_body`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ClothRobotArmSimulation | n/a | — | Constructed `ClothRobotArmSimulation`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:516-516`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Array shapes are only guaranteed by the constructor call site; no inner constructor checks that the time dimensions agree. Memory scales as 9 * n_joints * N_t plus the full state trajectory, which grows quickly for fine `dt_s`. Torques are stored in child-body frames, so comparing across joints requires rotating with the trajectory quaternions.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 36.
