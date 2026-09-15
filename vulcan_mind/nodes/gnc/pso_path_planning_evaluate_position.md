---
id: gnc.pso_path_planning_evaluate_position
label: evaluate_position
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: evaluate_position
  lines:
  - 368
  - 368
inputs:
- id: pos
  type: Any
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: cost_cutoff
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cost_cutoff` (default `Inf`).
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
  description: Return value of `evaluate_position`. Returns `rpo_normalized_path_cost_components(`.
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

# evaluate_position

## Purpose
Scores one flattened particle by converting it to a path and computing the normalised cost components.

## Design & Implementation
A closure that calls `rpo_position_to_path` with the captured start, goal and current waypoint count, then `rpo_normalized_path_cost_components` with the geometry, effective safety distance and an optional `cost_cutoff`. The cutoff lets the cost evaluator abandon a candidate early once it cannot beat the particle's personal best.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | Any | n/a | yes | Positional argument `pos`. |
| in | `cost_cutoff` | Any | n/a | no | Keyword argument `cost_cutoff` (default `Inf`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `evaluate_position`. Returns `rpo_normalized_path_cost_components(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:461-461`
- [[gnc.pso_path_planning_evaluate_swarm_bang|evaluate_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:409-409`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:368-368`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:268-268`

**Downstream**

- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:370-370`
- `callees` → [[gnc.pso_helpers_rpo_position_to_path|rpo_position_to_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:369-369`
<!-- vulcan:connections:end -->

## Limitations
Cost evaluation samples the whole path, so this is the dominant cost of the planner; `cost_cutoff` only helps when the evaluator honours it, which depends on the cost implementation.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 368.
