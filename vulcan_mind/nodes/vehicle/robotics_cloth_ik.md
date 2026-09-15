---
id: vehicle.robotics_cloth_ik
label: cloth_ik
kind: function
source:
  file: src/vehicle/robotics/robotics.jl
  symbol: cloth_ik
  lines:
  - 261
  - 261
inputs:
- id: model
  type: ClothArmModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: base_pose
  type: ClothArmBasePose
  units: n/a
  required: true
  description: Positional argument `base_pose`.
- id: target
  type: Any
  units: n/a
  required: true
  description: Positional argument `target`.
- id: q_seed
  type: Any
  units: n/a
  required: false
  description: Keyword argument `q_seed` (default `zeros(length(model.joints))`).
- id: max_iters
  type: Int
  units: n/a
  required: false
  description: Keyword argument `max_iters` (default `100`).
- id: position_tol_m
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `position_tol_m` (default `1.0e-4`).
- id: damping
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `damping` (default `1.0e-3`).
- id: step_limit_rad
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `step_limit_rad` (default `0.25`).
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
  type: Any
  units: n/a
  description: Return value of `cloth_ik`. Returns `q`.
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

# cloth_ik

## Purpose
Finds joint angles that place the end-effector at a target position, by damped least-squares iteration within joint limits.

## Theory & Math
Each iteration solves the Levenberg-Marquardt step

$$
\Delta q = \left( J^\top J + \lambda^2 I \right)^{-1} J^\top e,\qquad e = p_{\text{target}} - p_{\text{ee}}(q)
$$

with $J = \partial p_{\text{ee}} / \partial q$ and damping $\lambda$.

## Design & Implementation
Starts from `q_seed`, and for up to `max_iters` (100) iterations computes the position error, returning as soon as its norm is within `position_tol_m` (1e-4 m). Otherwise it forms the Jacobian and solves `(JᵀJ + λ²I) Δq = Jᵀ e` with `damping` 1e-3, scales `Δq` down to `step_limit_rad` (0.25) if longer, applies it and clamps every joint into its limits. After the loop it re-evaluates and throws `ErrorException` with the final error unless within ten times the tolerance, so a near-miss is accepted while a genuine failure is loud.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `target` | Any | n/a | yes | Positional argument `target`. |
| in | `q_seed` | Any | n/a | no | Keyword argument `q_seed` (default `zeros(length(model.joints))`). |
| in | `max_iters` | Int | n/a | no | Keyword argument `max_iters` (default `100`). |
| in | `position_tol_m` | Float64 | n/a | no | Keyword argument `position_tol_m` (default `1.0e-4`). |
| in | `damping` | Float64 | n/a | no | Keyword argument `damping` (default `1.0e-3`). |
| in | `step_limit_rad` | Float64 | n/a | no | Keyword argument `step_limit_rad` (default `0.25`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `cloth_ik`. Returns `q`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:188-188`
- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:99-99`
- [[grp.src_gnc_robotics|gnc/robotics/]] · `members_out` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:188-188`

**Downstream**

- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/vehicle/robotics/robotics.jl:274-274`
- `callees` → [[vehicle.robotics__ee_position_jacobian|_ee_position_jacobian]] · `callers` · call · `src/vehicle/robotics/robotics.jl:277-277`
- `callees` → [[vehicle.robotics__validate_joint_vector|_validate_joint_vector]] · `callers` · call · `src/vehicle/robotics/robotics.jl:271-271`
<!-- vulcan:connections:end -->

## Limitations
Only position is solved, never orientation, so the end-effector attitude is whatever the chain produces; the fixed damping does not adapt near singularities, and the clamp after each step can pin a joint at its limit and stall convergence.

## Provenance
Mapped from `src/vehicle/robotics/robotics.jl` line 261.
