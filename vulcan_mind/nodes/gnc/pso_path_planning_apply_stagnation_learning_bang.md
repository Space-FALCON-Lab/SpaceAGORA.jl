---
id: gnc.pso_path_planning_apply_stagnation_learning_bang
label: apply_stagnation_learning!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: apply_stagnation_learning!
  lines:
  - 440
  - 440
inputs:
- id: curr_cost
  type: Any
  units: n/a
  required: true
  description: Positional argument `curr_cost`.
- id: curr_obs
  type: Any
  units: n/a
  required: true
  description: Positional argument `curr_obs`.
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
  type: Nothing
  units: n/a
  description: Return value of `apply_stagnation_learning!`; mutates `curr_cost` in
    place. Returns `nothing`.
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

# apply_stagnation_learning!

## Purpose
Gives stalled particles a targeted nudge by copying one waypoint block at a time from the global best, accepting only moves that lower cost without worsening clearance.

## Design & Implementation
A closure that returns immediately unless learning is enabled with a positive threshold. It protects the elite fraction, caps trials per particle at `stagnation_learning_max_blocks` clamped to the waypoint count, and for each unprotected particle at or over the stagnation threshold tries waypoint blocks in a random order: copy the particle, overwrite one three-element block with `gbest`'s, evaluate with the current cost as cutoff, and accept if total improves by more than 1e-9 and `J_obs` does not worsen. Acceptance updates position, zeroes the block's velocity, and refreshes personal and global bests; the stagnation counter is then set through `rpo_pso_stagnation_count_after_learning`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `curr_cost` | Any | n/a | yes | Positional argument `curr_cost`. |
| in | `curr_obs` | Any | n/a | yes | Positional argument `curr_obs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `apply_stagnation_learning!`; mutates `curr_cost` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:440-440`

**Downstream**

- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:614-614`
- `callees` → [[gnc.planner_core_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:461-461`
- `callees` → [[gnc.pso_helpers_rpo_position_to_path|rpo_position_to_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:514-514`
- `callees` → [[gnc.pso_path_planning_cfg_for_current|cfg_for_current]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:611-611`
- `callees` → [[gnc.pso_path_planning_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:461-461`
- `callees` → [[gnc.pso_path_planning_evaluate_swarm_bang|evaluate_swarm!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:505-505`
- `callees` → [[gnc.pso_path_planning_iteration_timed_out|iteration_timed_out]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:533-533`
- `callees` → [[gnc.pso_path_planning_record_iteration_bang|record_iteration!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:527-527`
- `callees` → [[gnc.pso_path_planning_record_iteration_timeout_bang|record_iteration_timeout!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:528-528`
- `callees` → [[gnc.pso_path_planning_reexplore_waypoint_count|reexplore_waypoint_count]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:515-515`
- `callees` → [[gnc.pso_path_planning_reset_swarm_bang|reset_swarm!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:498-498`
- `callees` → [[gnc.pso_path_planning_rpo_pso_cull_swarm_bang|rpo_pso_cull_swarm!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:546-546`
- `callees` → [[gnc.pso_path_planning_rpo_pso_early_stopping_feasible|rpo_pso_early_stopping_feasible]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:575-575`
- `callees` → [[gnc.pso_path_planning_rpo_pso_iteration_weights|rpo_pso_iteration_weights]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:588-588`
- `callees` → [[gnc.pso_path_planning_rpo_pso_material_improvement|rpo_pso_material_improvement]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:569-569`
- `callees` → [[gnc.pso_path_planning_rpo_pso_protected_particle_mask|rpo_pso_protected_particle_mask]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:444-444`
- `callees` → [[gnc.pso_path_planning_rpo_pso_stagnation_count_after_learning|rpo_pso_stagnation_count_after_learning]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:493-493`
- `callees` → [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:613-613`
<!-- vulcan:connections:end -->

## Limitations
Each trial is a full path evaluation, so with many stagnant particles this phase can dominate an iteration and trip the runtime budget; the `candidate` buffer is shared across particles, which is safe only because the loop is serial.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 440.
