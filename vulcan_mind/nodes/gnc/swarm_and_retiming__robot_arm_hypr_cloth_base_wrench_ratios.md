---
id: gnc.swarm_and_retiming__robot_arm_hypr_cloth_base_wrench_ratios
label: _robot_arm_hypr_cloth_base_wrench_ratios
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_cloth_base_wrench_ratios
  lines:
  - 240
  - 240
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
  type: Tuple
  units: n/a
  description: Return value of `_robot_arm_hypr_cloth_base_wrench_ratios`. Returns
    `_robot_arm_hypr_rigid_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)`
    or `(force_ratio=force_ratio, torque_ratio=torque_ratio, node_ratio=node_ratio)`.
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

# _robot_arm_hypr_cloth_base_wrench_ratios

## Purpose
Estimates the peak base force and torque demand along a robot-arm trajectory by simulating the cloth-coupled arm dynamics and evaluating the coupled RHS at each sample.

## Design & Implementation
It first checks `isdefined(parentmodule(@__MODULE__), :ClothRobotArmDynamics)` and `:ClothMultibody`; if either is missing it delegates to the rigid estimator. Otherwise it resolves `simulate_cloth_robot_arm_plan`, `assign_coupled_cloth_robot_arm_rhs!` and `compliant_state_parts` via `getfield`, builds a `RobotArmPlan` with `_robot_arm_plan_from_q_reference(...; planner=:hypr)` targeting the FK end-effector position of `q_ref[:, end]`, and simulates it with `dt_s = cfg.retime_cloth_dt_s`, `duration_s = t_ref[end]`, `integrator = :implicit_midpoint`, and the four cloth stiffness/damping gains from `cfg`. For every sampled state it rebuilds a state NamedTuple through `_robot_arm_hypr_cloth_state_for_reaction`, allocates a zeroed `du`, and calls the RHS with `MVector{3}` accumulators `base_force` and `base_torque`. Ratios are `cfg.retime_base_wrench_model_gain * maximum(abs.(wrench)) / limit` for each finite limit; the sample is attributed to node `idx = clamp(searchsortedfirst(t_ref, t), 1, nt)` and also to `idx - 1`, taking the max.

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
| out | `result` | Tuple | n/a | — | Return value of `_robot_arm_hypr_cloth_base_wrench_ratios`. Returns `_robot_arm_hypr_rigid_base_wrench_ratios(model, base_pose, q_ref, t_ref, cfg)` or `(force_ratio=force_ratio, torque_ratio=torque_ratio, node_ratio=node_ratio)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.swarm_and_retiming__robot_arm_hypr_base_wrench_ratios|_robot_arm_hypr_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:341-341`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`

**Downstream**

- `callees` → [[gnc.planner_core__robot_arm_plan_from_q_reference|_robot_arm_plan_from_q_reference]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:258-258`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_state_for_reaction|_robot_arm_hypr_cloth_state_for_reaction]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:284-284`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_rigid_base_wrench_ratios|_robot_arm_hypr_rigid_base_wrench_ratios]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:249-249`
- `callees` → [[vehicle.cloth_fk|cloth_fk]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:257-257`
<!-- vulcan:connections:end -->

## Limitations
Uses the infinity norm of the wrench (`maximum(abs.(...))`), not the Euclidean norm, so ratios are direction dependent. Allocates a fresh `du` NamedTuple with several matrices per sample, which is expensive for long trajectories. Simulation failures are not caught. The spacecraft base state is hard-coded as a unit quaternion at the origin with `mass = 1.0`, so any base motion coupling is ignored.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 240.
