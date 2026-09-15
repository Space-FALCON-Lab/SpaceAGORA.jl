---
id: gnc.swarm_and_retiming__robot_arm_hypr_base_wrench_ratios
label: _robot_arm_hypr_base_wrench_ratios
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_base_wrench_ratios
  lines:
  - 333
  - 333
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
- id: q_ref
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `q_ref`.
- id: t_ref
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `t_ref`.
- id: cfg
  type: RobotArmHYPRConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_robot_arm_hypr_base_wrench_ratios`. Returns `_robot_arm_hypr_cloth_base_wrench_ratios(model,
    base_pose, q_ref, t_ref, cfg)` or `_robot_arm_hypr_rigid_base_wrench_ratios(model,
    base_pose, q_ref, t_ref, cfg)`.
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

# _robot_arm_hypr_base_wrench_ratios

## Purpose
Selects which reaction-load estimator to run for a candidate joint trajectory: the fast rigid-body finite-difference model or the full cloth-coupled multibody simulation.

## Design & Implementation
Signature `(model::ClothArmModel, base_pose::ClothArmBasePose, q_ref::Matrix{Float64}, t_ref::AbstractVector{<:Real}, cfg::RobotArmHYPRConfig)`. When `cfg.retime_cloth_physics_enable` is true it forwards to `_robot_arm_hypr_cloth_base_wrench_ratios`; otherwise to `_robot_arm_hypr_rigid_base_wrench_ratios`. Both return a NamedTuple `(force_ratio, torque_ratio, node_ratio)` where the ratios are demand divided by the configured base-wrench limits and `node_ratio` has one entry per column of `q_ref`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `q_ref` | Matrix{Float64} | n/a | yes | Positional argument `q_ref`. |
| in | `t_ref` | AbstractVector{<:Real} | n/a | yes | Positional argument `t_ref`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_base_wrench_ratios`. Returns `_robot_arm_hypr_cloth_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)` or `_robot_arm_hypr_rigid_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:340-340`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:205-205`
- [[gncz.swarm_and_retiming__robot_arm_hypr_retime_reference|_robot_arm_hypr_retime_reference]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:399-399`

**Downstream**

- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_base_wrench_ratios|_robot_arm_hypr_cloth_base_wrench_ratios]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:341-341`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_rigid_base_wrench_ratios|_robot_arm_hypr_rigid_base_wrench_ratios]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:343-343`
<!-- vulcan:connections:end -->

## Limitations
The choice is made purely on a boolean flag; the cloth path itself may still fall back to the rigid model if the optional `ClothRobotArmDynamics` and `ClothMultibody` modules are not loaded, and the caller is not told which model actually ran.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 333.
