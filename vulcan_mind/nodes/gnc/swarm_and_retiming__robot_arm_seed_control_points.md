---
id: gnc.swarm_and_retiming__robot_arm_seed_control_points
label: _robot_arm_seed_control_points
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_seed_control_points
  lines:
  - 85
  - 85
inputs:
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
- id: n_waypoints
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_waypoints`.
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
  description: Return value of `_robot_arm_seed_control_points`. Returns `points`.
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

# _robot_arm_seed_control_points

## Purpose
Generates the initial straight-line control-point matrix that seeds HYPR particles with a linear joint-space interpolation between start and goal.

## Design & Implementation
Converts `q_start` and `q_goal` to `Vector{Float64}`, allocates `points = zeros(n, n_waypoints + 2)`, and writes `q0` and `qf` into the first and last columns. Each interior column `j + 1` is `(1 - α) q0 + α qf` with `α = j / (n_waypoints + 1)`, giving `n_waypoints` evenly spaced waypoints strictly between the endpoints.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_start` | Any | n/a | yes | Positional argument `q_start`. |
| in | `q_goal` | Any | n/a | yes | Positional argument `q_goal`. |
| in | `n_waypoints` | Int | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_seed_control_points`. Returns `points`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:241-241`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Linear interpolation in joint space does not account for joint limits or obstacles; the seed may be infeasible and rely on the swarm to repair it. `q_start` and `q_goal` must be the same length. No wrap-around handling for revolute joints crossing ±π.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 85.
