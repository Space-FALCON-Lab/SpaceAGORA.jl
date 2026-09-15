---
id: gnc.swarm_and_retiming__robot_arm_flatten_internal_points
label: _robot_arm_flatten_internal_points
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_flatten_internal_points
  lines:
  - 100
  - 100
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
  description: Return value of `_robot_arm_flatten_internal_points`. Returns `flat`.
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

# _robot_arm_flatten_internal_points

## Purpose
Serialises the interior columns of a robot-arm control-point matrix into the flat `Vector{Float64}` layout that the HYPR particle swarm optimises over.

## Design & Implementation
Given `points` of size `n x (n_waypoints + 2)`, it infers `n = size(points, 1)` and `n_waypoints = size(points, 2) - 2`, allocates `flat = zeros(n * n_waypoints)`, and copies column `j + 1` for `j in 1:n_waypoints` into `flat[n*(j-1)+1 : n*j]`. Columns 1 and `end` (the fixed start and goal) are deliberately excluded. Layout is waypoint-major, joint-minor, matching `_robot_arm_control_points`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_flatten_internal_points`. Returns `flat`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:245-245`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A matrix with fewer than two columns yields a negative `n_waypoints` and `zeros` throws `ArgumentError`. The function assumes `points` is a dense column-indexable matrix; it allocates once and does not accept a preallocated output buffer.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 100.
