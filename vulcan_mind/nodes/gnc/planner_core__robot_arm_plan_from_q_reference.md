---
id: gnc.planner_core__robot_arm_plan_from_q_reference
label: _robot_arm_plan_from_q_reference
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/planner_core.jl
  symbol: _robot_arm_plan_from_q_reference
  lines:
  - 121
  - 121
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
- id: q_start
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_start`.
- id: q_goal
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_goal`.
- id: target
  type: Any
  units: n/a
  required: true
  description: Positional argument `target`.
- id: t_ref
  type: Any
  units: n/a
  required: true
  description: Positional argument `t_ref`.
- id: q_ref
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_ref`.
- id: planner
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `planner` (default `:hypr`).
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
  type: RobotArmPlan
  units: n/a
  description: Return value of `_robot_arm_plan_from_q_reference`. Returns `RobotArmPlan(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _robot_arm_plan_from_q_reference

## Purpose
Packs a joint-space reference `q_ref` (joints by samples, rad) and its timestamps `t_ref` (s) into a `RobotArmPlan`, filling in numerically differentiated joint velocities and accelerations, the end-effector trajectory from forward kinematics, and the final position error against the Cartesian `target`.

## Theory & Math
For non-uniform spacing with $h_- = t_k - t_{k-1}$ and $h_+ = t_{k+1} - t_k$:
$$\dot q_k \approx \frac{q_{k+1} - q_{k-1}}{h_- + h_+}, \qquad \ddot q_k \approx \frac{2}{h_- + h_+}\left(\frac{q_{k+1} - q_k}{h_+} - \frac{q_k - q_{k-1}}{h_-}\right)$$

## Design & Implementation
`dq_ref` uses one-sided differences at the first and last sample and a non-uniform central difference `(q[k+1] - q[k-1]) / (dt_prev + dt_next)` in the interior. `ddq_ref` uses the non-uniform second difference `2*((q[k+1]-q[k])/dt_next - (q[k]-q[k-1])/dt_prev)/(dt_prev+dt_next)` and copies the neighbouring interior value into the two end samples. `ee_ref` is built column by column from `cloth_fk(model, base_pose, q_ref[:, k]).end_effector_position`. `final_error = norm(ee_ref[:, end] - target)`. The constructor receives `collect(Float64.(t_ref))`, `Float64` copies of `q_start` and `q_goal`, the target as `SVector{3,Float64}`, and the `planner` symbol (default `:hypr`).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `q_start` | Any | n/a | yes | Positional argument `q_start`. |
| in | `q_goal` | Any | n/a | yes | Positional argument `q_goal`. |
| in | `target` | Any | n/a | yes | Positional argument `target`. |
| in | `t_ref` | Any | n/a | yes | Positional argument `t_ref`. |
| in | `q_ref` | Any | n/a | yes | Positional argument `q_ref`. |
| in | `planner` | Symbol | n/a | no | Keyword argument `planner` (default `:hypr`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmPlan | n/a | — | Return value of `_robot_arm_plan_from_q_reference`. Returns `RobotArmPlan(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:339-339`
- [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_base_wrench_ratios|_robot_arm_hypr_cloth_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:258-258`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:204-204`

**Downstream**

- `callees` → [[gnc.robot_arm_planning_robotarmplan|RobotArmPlan]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:157-157`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:153-153`
<!-- vulcan:connections:end -->

## Limitations
If `nt == 1`, velocities and accelerations are all zero; if `nt == 2`, `ddq_ref` stays zero. Repeated timestamps produce division by zero and `Inf`/`NaN` derivatives with no check. The central difference for velocity is only first-order accurate on non-uniform grids. Forward kinematics is evaluated once per sample, so cost scales linearly with `nt`. `q_ref` is stored by reference, not copied, so later mutation by the caller changes the plan.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/planner_core.jl` line 121.
