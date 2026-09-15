---
id: gnc.swarm_and_retiming_robot_arm_sample_hypr_path
label: robot_arm_sample_hypr_path
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: robot_arm_sample_hypr_path
  lines:
  - 112
  - 112
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: n_samples
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_samples`.
- id: curve_type
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `curve_type` (default `:bezier`).
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
  description: Return value of `robot_arm_sample_hypr_path`. Returns `hypr_sample_count_path(points,
    n_samples; curve_type=curve_type)`.
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

# robot_arm_sample_hypr_path

## Purpose
Public entry that samples a robot-arm control-point matrix into a dense joint-space path, either as a Bezier curve or as a straight polyline, for cost evaluation and execution.

## Design & Implementation
Signature `robot_arm_sample_hypr_path(points, n_samples::Int; curve_type::Symbol = :bezier)`. It forwards to the shared `hypr_sample_count_path(points, n_samples; curve_type)`, which returns a joints x `n_samples` matrix. `points` is the `n x (n_waypoints + 2)` matrix from `_robot_arm_control_points` or `_robot_arm_seed_control_points`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `n_samples` | Int | n/a | yes | Positional argument `n_samples`. |
| in | `curve_type` | Symbol | n/a | no | Keyword argument `curve_type` (default `:bezier`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `robot_arm_sample_hypr_path`. Returns `hypr_sample_count_path(points, n_samples; curve_type=curve_type)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:335-335`
- [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:10-10`
- [[gnc.rrt_warmstart__robot_arm_rrt_path_score|_robot_arm_rrt_path_score]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:176-176`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:200-200`

**Downstream**

- `callees` → [[gncz.hypr_utils_hypr_sample_count_path|hypr_sample_count_path]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:113-113`
<!-- vulcan:connections:end -->

## Limitations
Accepted values of `curve_type` are defined by `hypr_sample_count_path`; an unsupported symbol errors there rather than here. A Bezier curve of high degree does not pass through its interior control points, so waypoint constraints are only approximately honoured.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 112.
