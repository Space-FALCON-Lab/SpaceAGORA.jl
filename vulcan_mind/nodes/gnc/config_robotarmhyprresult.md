---
id: gnc.config_robotarmhyprresult
label: RobotArmHYPRResult
kind: struct
source:
  file: src/gnc/robotics/robot_arm_hypr/config.jl
  symbol: RobotArmHYPRResult
  lines:
  - 79
  - 79
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Field `plan`.
- id: control_points
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `control_points`.
- id: sampled_path
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `sampled_path`.
- id: cost
  type: Float64
  units: n/a
  required: true
  description: Field `cost`.
- id: components
  type: NamedTuple
  units: n/a
  required: true
  description: Field `components`.
- id: cost_history
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `cost_history`.
- id: config
  type: RobotArmHYPRConfig
  units: n/a
  required: true
  description: Field `config`.
- id: early_stopped
  type: Bool
  units: n/a
  required: true
  description: Field `early_stopped`.
- id: early_stop_iter
  type: Int
  units: n/a
  required: true
  description: Field `early_stop_iter`.
- id: cull_replacements
  type: Int
  units: n/a
  required: true
  description: Field `cull_replacements`.
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
  type: RobotArmHYPRResult
  units: n/a
  description: Constructed `RobotArmHYPRResult`.
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

# RobotArmHYPRResult

## Purpose
`RobotArmHYPRResult` is the return record of a robot-arm HYPR (hybrid PSO) planning run. It bundles the executable plan with everything needed to audit or plot how the optimiser arrived at it: the optimised control points, the densely sampled joint path, the final scalar cost and its breakdown, the per-iteration cost history, the exact configuration used, and the early-stopping and culling diagnostics.

## Design & Implementation
An immutable struct with ten fields: `plan::RobotArmPlan`, `control_points::Matrix{Float64}` and `sampled_path::Matrix{Float64}` (joint angles in radians), `cost::Float64`, `components::NamedTuple` holding the weighted length, smoothness and obstacle terms, `cost_history::Vector{Float64}` with one entry per completed iteration, `config::RobotArmHYPRConfig` echoed back verbatim, and the diagnostics `early_stopped::Bool`, `early_stop_iter::Int` and `cull_replacements::Int`. Echoing the configuration makes a result self-describing, so a saved result can be replayed or compared without the caller tracking which settings produced it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Field `plan`. |
| in | `control_points` | Matrix{Float64} | n/a | yes | Field `control_points`. |
| in | `sampled_path` | Matrix{Float64} | n/a | yes | Field `sampled_path`. |
| in | `cost` | Float64 | n/a | yes | Field `cost`. |
| in | `components` | NamedTuple | n/a | yes | Field `components`. |
| in | `cost_history` | Vector{Float64} | n/a | yes | Field `cost_history`. |
| in | `config` | RobotArmHYPRConfig | n/a | yes | Field `config`. |
| in | `early_stopped` | Bool | n/a | yes | Field `early_stopped`. |
| in | `early_stop_iter` | Int | n/a | yes | Field `early_stop_iter`. |
| in | `cull_replacements` | Int | n/a | yes | Field `cull_replacements`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmHYPRResult | n/a | — | Constructed `RobotArmHYPRResult`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:352-352`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:218-218`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`components` is an untyped `NamedTuple`, so its field names and arity are a convention held only by the producing code — consumers that index it by name break silently if the cost decomposition changes. The struct is immutable but its `Matrix` and `Vector` fields are not, so a caller can mutate `sampled_path` or `cost_history` in place and desynchronise them from `cost`. Nothing in the record states whether the returned path is collision-free; feasibility must be re-derived from `components`.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/config.jl` line 79.
