---
id: simulation.dynamics_rhs__count_non_batchable_effectors
label: _count_non_batchable_effectors
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _count_non_batchable_effectors
  lines:
  - 674
  - 674
inputs:
- id: effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `effectors`.
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
  type: Int
  units: n/a
  description: Return value of `_count_non_batchable_effectors`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _count_non_batchable_effectors

## Purpose
Counts effectors without a vectorised batch kernel, used by the planner to judge whether the flat queue has enough generic work to be worthwhile.

## Design & Implementation
Loops the tuple and counts effectors failing `_batchable_effector`. Declared `@inline` with an `::Int` return. Differs from the flat-queue-only count in that harmonics is counted as non-batchable here even though it has its own prepass.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_count_non_batchable_effectors`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1047-1047`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
None.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 674.
