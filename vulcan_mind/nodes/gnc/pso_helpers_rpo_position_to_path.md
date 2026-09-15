---
id: gnc.pso_helpers_rpo_position_to_path
label: rpo_position_to_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_helpers.jl
  symbol: rpo_position_to_path
  lines:
  - 25
  - 25
inputs:
- id: position
  type: Any
  units: n/a
  required: true
  description: Positional argument `position`.
- id: start_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `start_rtn`.
- id: goal_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal_rtn`.
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
  description: Return value of `rpo_position_to_path`. Returns `points`.
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

# rpo_position_to_path

## Purpose
Decodes a flattened particle swarm design vector into the full three-by-N waypoint matrix of a candidate relative-motion path, with the fixed start and goal reattached at the ends.

## Design & Implementation
Allocates `points = zeros(3, n_waypoints + 2)`, writes `SVector{3, Float64}(start_rtn)` into column one and `SVector{3, Float64}(goal_rtn)` into the last column, then for each internal waypoint `j` in `1:n_waypoints` copies `position[(3 * (j - 1) + 1):(3 * j)]` into column `j + 1`. The encoding is therefore contiguous triples in RTN metres, waypoint-major, and the boundary conditions are never part of the optimised vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `position` | Any | n/a | yes | Positional argument `position`. |
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `n_waypoints` | Int | n/a | yes | Positional argument `n_waypoints`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_position_to_path`. Returns `points`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:514-514`
- [[gnc.pso_path_planning_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:369-369`
- [[gnc.pso_path_planning_rpo_pso_cull_swarm_bang|rpo_pso_cull_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:77-77`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:369-369`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No check couples `length(position)` to `3 * n_waypoints`; a short vector raises a bounds error mid-loop after partially filling the matrix, and a longer one silently ignores its tail. A fresh dense matrix is allocated on every call, which is significant because the cost function evaluates it once per particle per iteration. Waypoint ordering is taken as given, so the decoded path may self-intersect or backtrack.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_helpers.jl` line 25.
