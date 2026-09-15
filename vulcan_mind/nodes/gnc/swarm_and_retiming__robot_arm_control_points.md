---
id: gnc.swarm_and_retiming__robot_arm_control_points
label: _robot_arm_control_points
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_control_points
  lines:
  - 70
  - 70
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
- id: pos
  type: Any
  units: n/a
  required: true
  description: Positional argument `pos`.
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
  description: Return value of `_robot_arm_control_points`. Returns `points`.
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

# _robot_arm_control_points

## Purpose
Rebuilds the full joint-space control-point matrix for a robot-arm HYPR particle by bracketing its flattened internal waypoints with the fixed start and goal configurations.

## Design & Implementation
Takes `q_start`, `q_goal` (any iterable of joint angles, converted via `Float64.(collect(...))`), the flat particle vector `pos`, and `n_waypoints::Int`. It allocates `points = zeros(n, n_waypoints + 2)` where `n = length(q0)`, writes `q0` into column 1 and `qf` into the last column, then for each waypoint `j` copies the slice `pos[n*(j-1)+1 : n*j]` into column `j + 1` using a `@view` under `@inbounds`. This is the exact inverse of `_robot_arm_flatten_internal_points`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_start` | Any | n/a | yes | Positional argument `q_start`. |
| in | `q_goal` | Any | n/a | yes | Positional argument `q_goal`. |
| in | `pos` | Any | n/a | yes | Positional argument `pos`. |
| in | `n_waypoints` | Int | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_control_points`. Returns `points`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:269-269`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:269-269`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No check that `length(pos) == n * n_waypoints`; a shorter `pos` reads out of bounds silently under `@inbounds`, and a longer one is truncated. `q_start` and `q_goal` are assumed to have equal length; a mismatch throws a `DimensionMismatch` from the broadcast on the last column. Allocates a fresh matrix per call, which matters inside the PSO cost loop.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 70.
