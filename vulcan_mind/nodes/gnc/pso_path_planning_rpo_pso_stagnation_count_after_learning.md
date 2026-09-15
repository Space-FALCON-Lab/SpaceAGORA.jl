---
id: gnc.pso_path_planning_rpo_pso_stagnation_count_after_learning
label: rpo_pso_stagnation_count_after_learning
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_stagnation_count_after_learning
  lines:
  - 179
  - 179
inputs:
- id: count
  type: Integer
  units: n/a
  required: true
  description: Positional argument `count`.
- id: accepted
  type: Bool
  units: n/a
  required: true
  description: Positional argument `accepted`.
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
  description: 'Return value of `rpo_pso_stagnation_count_after_learning`. Returns
    `accepted ? 0 : Int(count)`.'
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

# rpo_pso_stagnation_count_after_learning

## Purpose
Defines the stagnation counter's behaviour after a learning attempt: only an accepted move clears it.

## Design & Implementation
An `@inline` conditional returning zero when `accepted` and the unchanged count otherwise. Factoring the rule into a named function makes the semantics testable in isolation and documents that a rejected attempt does not count as progress.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `count` | Integer | n/a | yes | Positional argument `count`. |
| in | `accepted` | Bool | n/a | yes | Positional argument `accepted`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_stagnation_count_after_learning`. Returns `accepted ? 0 : Int(count)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:493-493`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:493-493`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A particle whose attempts are repeatedly rejected keeps its high count and is re-tried every iteration, spending evaluations without a back-off.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 179.
