---
id: gnc.robot_arm_planning_robotarmplanning
label: RobotArmPlanning
kind: module
source:
  file: src/gnc/robotics/robot_arm_planning.jl
  symbol: RobotArmPlanning
  lines:
  - 2
  - 2
inputs:
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
  description: Value produced by this symbol.
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

# RobotArmPlanning

## Purpose
Module providing joint-space motion planning for the cloth-multibody robot arm: a quintic time-scaled interpolation planner between an initial joint vector and an IK-solved goal, plus an included HYPR sampling-based planner with sphere obstacles.

## Design & Implementation
`RobotArmPlanning` imports `LinearAlgebra`, `Random`, `StaticArrays`, the sibling `Robotics` module (for `ClothArmModel`, `ClothArmBasePose`, `cloth_ik`, `cloth_fk`) and `HYPRUtils`. It defines `RobotArmPlannerConfig`, `RobotArmPlan`, the private quintic blend functions `_quintic_scalar`, `_quintic_scalar_dot`, `_quintic_scalar_ddot`, the grid builder `_reference_times`, then `include`s `robot_arm_hypr.jl` before defining `plan_robot_arm_motion` (dispatching on `planner=:cloth_quintic` or `:hypr`, throwing `ArgumentError` for anything else) and `robot_arm_plan_sample`. Exports cover the config/plan types, both planners, the HYPR config/result/obstacle types and the HYPR sampling and cost helpers.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_planning.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The quintic planner is purely kinematic: it ignores joint limits, obstacles and dynamics, and relies entirely on `cloth_ik` converging (a failed IK still yields a plan whose `final_error_m` is simply large). Obstacles and `hypr_config` are accepted by `plan_robot_arm_motion` but silently ignored unless `planner == :hypr`. The include order means `robot_arm_hypr.jl` must not reference `plan_robot_arm_motion` at load time.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_planning.jl` line 2.
