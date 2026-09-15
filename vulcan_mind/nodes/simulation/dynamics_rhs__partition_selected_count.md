---
id: simulation.dynamics_rhs__partition_selected_count
label: _partition_selected_count
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _partition_selected_count
  lines:
  - 42
  - 42
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: partition
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `partition`.
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
  description: Return value of `_partition_selected_count`.
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

# _partition_selected_count

## Purpose
Counts how many effectors belong to a given solver partition, letting the partitioned accumulators return early when the partition is empty and decide whether inner effector threading could pay off.

## Design & Implementation
Loops the effector tuple under `@inbounds` summing one for each effector where `_effector_in_partition` holds. Declared `@inline` with an `::Int` return. The threaded path in `_accumulate_dynamic_effectors_partitioned!` requires this count to exceed one.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `partition` | Symbol | n/a | yes | Positional argument `partition`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_partition_selected_count`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1004-1004`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:125-125`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Recomputed on every partitioned RHS call although the effector tuple is fixed for a run; caching it at setup would remove a small constant cost.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 42.
