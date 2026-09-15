---
id: dynamics.cloth_multibody_compliant_joint_loads
label: compliant_joint_loads
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: compliant_joint_loads
  lines:
  - 528
  - 528
inputs:
- id: model
  type: CompliantMultibodyModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `x`.
- id: joint_rest_quaternions
  type: Union{Nothing, AbstractVector}
  units: n/a
  required: false
  description: Keyword argument `joint_rest_quaternions` (default `nothing`).
- id: joint_actuators
  type: AbstractVector{CompliantJointActuator}
  units: n/a
  required: false
  description: Keyword argument `joint_actuators` (default `CompliantJointActuator[]`).
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
  type: Vector{CompliantJointLoad}
  units: n/a
  description: Return value of `compliant_joint_loads`.
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

# compliant_joint_loads

## Purpose
Computes the spring-damper and actuator loads at every joint for a given state, both to drive the dynamics and to expose diagnostics.

## Theory & Math
$$
\vec{F}_{p} = -K_t (\vec{p}_p - \vec{p}_c) - C_t (\vec{v}_p - \vec{v}_c),\qquad \vec{\tau}_{c} = K_r \vec{\phi} - C_r (\vec{\omega}_c - \vec{\omega}_p)
$$

## Design & Implementation
Validates the state length. For each joint it gathers both sides' kinematics, forms the attachment-point separation and relative velocity, and applies the translational law `-K Δr - C Δv`. The rotational law compares the child orientation to `parent_q ⊗ rest` via `_axis_angle_error` and applies `K_rot φ - C_rot Δω`; actuators bound to the joint add their torques. Parent-side loads are the negatives, and child-body-frame torques are provided. Returns a vector of `CompliantJointLoad`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | CompliantMultibodyModel | n/a | yes | Positional argument `model`. |
| in | `x` | AbstractVector | n/a | yes | Positional argument `x`. |
| in | `joint_rest_quaternions` | Union{Nothing, AbstractVector} | n/a | no | Keyword argument `joint_rest_quaternions` (default `nothing`). |
| in | `joint_actuators` | AbstractVector{CompliantJointActuator} | n/a | no | Keyword argument `joint_actuators` (default `CompliantJointActuator[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{CompliantJointLoad} | n/a | — | Return value of `compliant_joint_loads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:504-504`
- [[dynx.multibody_cloth_compliant_multibody_dynamics|compliant_multibody_dynamics]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:595-595`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__actuator_torque_child_world|_actuator_torque_child_world]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:555-555`
- `callees` → [[dynamics.cloth_multibody__axis_angle_error|_axis_angle_error]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:549-549`
- `callees` → [[dynamics.cloth_multibody__joint_rest|_joint_rest]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:546-546`
- `callees` → [[dynamics.cloth_multibody__parent_kinematics|_parent_kinematics]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:538-538`
- `callees` → [[dynamics.cloth_multibody__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:548-548`
- `callees` → [[dynamics.cloth_multibody__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:547-547`
- `callees` → [[dynamics.cloth_multibody_compliantjointload|CompliantJointLoad]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:562-562`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__axis_angle_error|_axis_angle_error]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:549-549`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_conj|_quat_conj]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:548-548`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:547-547`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:562-562`
- `callees` → [[vehicle.robotics__quat_mul|_quat_mul]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:547-547`
<!-- vulcan:connections:end -->

## Limitations
Every `_quat_mul` normalises, and every joint allocates a load struct; for a large grid this is the dominant cost of the derivative.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 528.
