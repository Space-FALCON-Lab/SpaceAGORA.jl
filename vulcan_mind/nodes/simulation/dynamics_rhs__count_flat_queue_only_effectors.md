---
id: simulation.dynamics_rhs__count_flat_queue_only_effectors
label: _count_flat_queue_only_effectors
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _count_flat_queue_only_effectors
  lines:
  - 695
  - 695
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
  description: Return value of `_count_flat_queue_only_effectors`.
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

# _count_flat_queue_only_effectors

## Purpose
Counts effectors that must go through the generic dynamic flat queue because they have neither a vectorised batch kernel nor the harmonics prepass, which decides whether per-worker partial buffers need zeroing.

## Design & Implementation
Loops the effector tuple and increments for each that fails both `_batchable_effector` and `_harmonics_prepass_effector`. Declared `@inline` with an `::Int` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_count_flat_queue_only_effectors`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1002-1002`
- [[simulation.dynamics_rhs__prepare_rhs_flat_work_items_bang|_prepare_rhs_flat_work_items!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:429-429`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
None beyond depending on the two predicates staying in sync with the batch dispatch.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 695.
