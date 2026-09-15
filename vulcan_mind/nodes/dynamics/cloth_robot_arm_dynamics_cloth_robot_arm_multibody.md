---
id: dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody
label: cloth_robot_arm_multibody
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl
  symbol: cloth_robot_arm_multibody
  lines:
  - 225
  - 225
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: k_translation_n_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `k_translation_n_m` (default `5.0e3`).
- id: c_translation_n_s_m
  type: Any
  units: n/a
  required: false
  description: Keyword argument `c_translation_n_s_m` (default `30.0`).
- id: k_rotation_n_m_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `k_rotation_n_m_rad` (default `15.0`).
- id: c_rotation_n_m_s_rad
  type: Any
  units: n/a
  required: false
  description: Keyword argument `c_rotation_n_m_s_rad` (default `0.5`).
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
  type: CompliantMultibodyModel
  units: n/a
  description: Return value of `cloth_robot_arm_multibody`. Returns `CompliantMultibodyModel(bodies,
    joints, plan.base_pose.position, plan.base_pose.`.
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

# cloth_robot_arm_multibody

## Purpose
Builds the `CompliantMultibodyModel` for a robot-arm plan: one cylindrical `CompliantBody` per link and one spring-damper `CompliantJoint` per link connecting it to its parent (the base for link 1).

## Design & Implementation
Evaluates `cloth_fk` at `plan.q_ref[:,1]` for parent quaternions, constructs bodies with `_link_inertia`, and expands the four keyword compliances (`k_translation_n_m=5e3`, `c_translation_n_s_m=30`, `k_rotation_n_m_rad=15`, `c_rotation_n_m_s_rad=0.5`) via `_joint_compliance_matrices`. Rest quaternions come from `cloth_robot_arm_rest_quaternions(plan, 0.0)`. For joint i, the parent attachment point is `mount_offset_body` (i == 1) or `vector_parent - com_offset_parent` of link i-1, and the child point is `-com_offset_parent`. Joints are named `:cloth_arm_joint_i`. The model is returned with `plan.base_pose.position` and `.quaternion` as the base.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `k_translation_n_m` | Any | n/a | no | Keyword argument `k_translation_n_m` (default `5.0e3`). |
| in | `c_translation_n_s_m` | Any | n/a | no | Keyword argument `c_translation_n_s_m` (default `30.0`). |
| in | `k_rotation_n_m_rad` | Any | n/a | no | Keyword argument `k_rotation_n_m_rad` (default `15.0`). |
| in | `c_rotation_n_m_s_rad` | Any | n/a | no | Keyword argument `c_rotation_n_m_s_rad` (default `0.5`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantMultibodyModel | n/a | — | Return value of `cloth_robot_arm_multibody`. Returns `CompliantMultibodyModel(bodies, joints, plan.base_pose.position, plan.base_pose.`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:451-451`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody_compliantbody|CompliantBody]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:235-235`
- `callees` → [[dynamics.cloth_multibody_compliantjoint|CompliantJoint]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:257-257`
- `callees` → [[dynamics.cloth_multibody_compliantmultibodymodel|CompliantMultibodyModel]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:270-270`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__joint_compliance_matrices|_joint_compliance_matrices]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:238-238`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__link_inertia|_link_inertia]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:235-235`
- `callees` → [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_rest_quaternions|cloth_robot_arm_rest_quaternions]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:242-242`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:257-257`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:233-233`
<!-- vulcan:connections:end -->

## Limitations
Rest quaternions are frozen at t = 0 inside the model even though the steppers accept per-step `joint_rest_quaternions`; the model's own values are only meaningful for the initial configuration. Assumes a strictly serial chain (parent of link i is link i-1). `pose0` is computed but only its quaternions are used.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl` line 225.
