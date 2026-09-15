---
id: gnc.planner_comparison_rpo_flatten_planner_results
label: rpo_flatten_planner_results
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_flatten_planner_results
  lines:
  - 649
  - 649
inputs:
- id: batch
  type: Any
  units: n/a
  required: true
  description: Positional argument `batch`.
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
  description: Return value of `rpo_flatten_planner_results`. Returns `rows`.
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

# rpo_flatten_planner_results

## Purpose
Flattens the per-planner result dictionary from a comparison batch into a single vector of NamedTuple rows, in planner order, ready for CSV export and metric aggregation.

## Design & Implementation
`rpo_flatten_planner_results(batch)` iterates `batch.planner_types` and `append!`s `batch.results_by_planner[planner]` to a `NamedTuple[]` accumulator, returning it. Planner order follows `planner_types`, which may have been reordered to place `:hypr` first when runtime matching is enabled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `batch` | Any | n/a | yes | Positional argument `batch`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_flatten_planner_results`. Returns `rows`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_write_planner_comparison_outputs|rpo_write_planner_comparison_outputs]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1156-1156`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The container is `Vector{NamedTuple}` (abstract element type), so downstream field access is dynamically dispatched. It assumes every planner in `planner_types` has an entry in `results_by_planner`, which holds for batches produced by `rpo_run_planner_comparison_batch` but raises `KeyError` for hand-built batches.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 649.
