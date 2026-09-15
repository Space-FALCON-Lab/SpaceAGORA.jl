---
id: gnc.robot_arm_planning_robotarmplan
label: RobotArmPlan
kind: struct
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: RobotArmPlan
  lines:
  - 25
  - 25
inputs:
- id: model
  type: ClothArmModel
  units: n/a
  required: true
  description: Field `model`.
- id: base_pose
  type: ClothArmBasePose
  units: n/a
  required: true
  description: Field `base_pose`.
- id: t_ref_s
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `t_ref_s`.
- id: q_ref
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `q_ref`.
- id: dq_ref
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `dq_ref`.
- id: ddq_ref
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `ddq_ref`.
- id: ee_ref
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `ee_ref`.
- id: q_start
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `q_start`.
- id: q_goal
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `q_goal`.
- id: target
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `target`.
- id: final_error_m
  type: Float64
  units: n/a
  required: true
  description: Field `final_error_m`.
- id: planner
  type: Symbol
  units: n/a
  required: true
  description: Field `planner`.
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
  description: Constructed `RobotArmPlan`.
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

# RobotArmPlan

## Purpose
Immutable container holding a planned robot-arm motion: the arm model, base pose, time grid, joint position/velocity/acceleration references, end-effector trajectory, boundary joint vectors, target point, achieved final error and which planner produced it.

## Design & Implementation
Fields are `model::ClothArmModel`, `base_pose::ClothArmBasePose`, `t_ref_s::Vector{Float64}` (s), `q_ref`, `dq_ref`, `ddq_ref::Matrix{Float64}` each sized `n_joints x n_times` (rad, rad/s, rad/s^2), `ee_ref::Matrix{Float64}` sized `3 x n_times` (m), `q_start` and `q_goal::Vector{Float64}` (rad), `target::SVector{3,Float64}` (m), `final_error_m::Float64` (norm of the last `ee_ref` column minus `target`) and `planner::Symbol` (`:cloth_quintic` or `:hypr`). It is constructed positionally by `plan_robot_arm_motion` and the HYPR planner and consumed by `robot_arm_plan_sample`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Field `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Field `base_pose`. |
| in | `t_ref_s` | Vector{Float64} | n/a | yes | Field `t_ref_s`. |
| in | `q_ref` | Matrix{Float64} | n/a | yes | Field `q_ref`. |
| in | `dq_ref` | Matrix{Float64} | n/a | yes | Field `dq_ref`. |
| in | `ddq_ref` | Matrix{Float64} | n/a | yes | Field `ddq_ref`. |
| in | `ee_ref` | Matrix{Float64} | n/a | yes | Field `ee_ref`. |
| in | `q_start` | Vector{Float64} | n/a | yes | Field `q_start`. |
| in | `q_goal` | Vector{Float64} | n/a | yes | Field `q_goal`. |
| in | `target` | SVector{3, Float64} | n/a | yes | Field `target`. |
| in | `final_error_m` | Float64 | n/a | yes | Field `final_error_m`. |
| in | `planner` | Symbol | n/a | yes | Field `planner`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmPlan | n/a | — | Constructed `RobotArmPlan`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core__robot_arm_plan_from_q_reference|_robot_arm_plan_from_q_reference]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:157-157`
- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:128-128`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct performs no shape validation, so mismatched column counts between `t_ref_s` and the reference matrices only surface as bounds errors during sampling. Matrices are stored column-major per time sample and are mutable despite the outer struct being immutable, so a caller can corrupt a plan in place. `final_error_m` reflects IK convergence only for the goal sample.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 25.
