---
id: gnc.pso_path_planning_record_iteration_bang
label: record_iteration!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: record_iteration!
  lines:
  - 265
  - 265
inputs:
- id: iter
  type: Any
  units: n/a
  required: true
  description: Positional argument `iter`.
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
  description: Return value of `record_iteration!`; mutates `iter` in place. Returns
    `nothing`.
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

# record_iteration!

## Purpose
Appends the current global best cost to the history and notifies the optional iteration callback.

## Design & Implementation
A closure pushing `gbest_cost` onto `cost_history` and, if `iteration_callback` is not `nothing`, invoking it with the iteration index, cost and best components. It is called once per completed iteration and once on each timeout path before breaking, so the history length equals the number of iterations that ran.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `iter` | Any | n/a | yes | Positional argument `iter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `record_iteration!`; mutates `iter` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:527-527`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:265-265`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:266-266`
<!-- vulcan:connections:end -->

## Limitations
The callback runs synchronously on the planner's thread, so a slow callback counts against the iteration runtime budget.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 265.
