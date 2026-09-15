---
id: gnc.pso_path_planning_reexplore_waypoint_count
label: reexplore_waypoint_count
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: reexplore_waypoint_count
  lines:
  - 431
  - 431
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
  description: Return value of `reexplore_waypoint_count`. Returns `target`.
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

# reexplore_waypoint_count

## Purpose
Decides how many waypoints to grow to when the best path is still in violation and reexploration triggers.

## Design & Implementation
A closure taking the largest of the current count plus one, the count scaled by `reexplore_waypoint_scale` rounded up, and the count plus `reexplore_waypoint_increment`, then capping at `reexplore_max_waypoints` when that is positive. Growing by at least one guarantees progress each trigger.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `reexplore_waypoint_count`. Returns `target`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:515-515`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:431-431`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Once the cap is reached the count stops growing but the search margin can still expand, so repeated triggers thereafter only widen the box.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 431.
