---
id: gnc.pso_path_planning_evaluate_swarm_bang
label: evaluate_swarm!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: evaluate_swarm!
  lines:
  - 402
  - 402
inputs:
- id: iter_start_ns
  type: Any
  units: n/a
  required: false
  description: Keyword argument `iter_start_ns` (default `nothing`).
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
  description: Return value of `evaluate_swarm!`; mutates `iter_start_ns` in place.
    Returns `curr_cost, curr_obs, iteration_timed_out(iter_start_ns)` or `curr_cost,
    curr_obs, false`.
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

# evaluate_swarm!

## Purpose
Evaluates every particle for the current iteration and updates the best-state bookkeeping, honouring the runtime budget when one is set.

## Design & Implementation
A closure. When `iter_start_ns` is given and a finite positive iteration limit exists, it evaluates particles serially, checking `iteration_timed_out` before each so it can return a partial result with the timed-out flag. Otherwise it evaluates in parallel with `@threads`, storing components per particle, then applies `update_swarm_best_from_components!` serially. It returns the current cost vector, obstacle-cost vector and the timeout flag.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `iter_start_ns` | Any | n/a | no | Keyword argument `iter_start_ns` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `evaluate_swarm!`; mutates `iter_start_ns` in place. Returns `curr_cost, curr_obs, iteration_timed_out(iter_start_ns)` or `curr_cost, curr_obs, false`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:505-505`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:402-402`

**Downstream**

- `callees` → [[gnc.planner_core_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:409-409`
- `callees` → [[gnc.pso_path_planning_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:409-409`
- `callees` → [[gnc.pso_path_planning_iteration_timed_out|iteration_timed_out]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:408-408`
- `callees` → [[gnc.pso_path_planning_update_swarm_best_from_components_bang|update_swarm_best_from_components!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:413-413`
<!-- vulcan:connections:end -->

## Limitations
Under a runtime limit the evaluation is single-threaded, so enabling the budget also disables the main source of parallelism; a partial evaluation leaves the unevaluated particles' entries at `Inf`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 402.
