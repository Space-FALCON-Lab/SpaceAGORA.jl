---
id: simulation.dynamics_rhs__partition_needs_state_sample
label: _partition_needs_state_sample
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _partition_needs_state_sample
  lines:
  - 35
  - 35
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
  type: Bool
  units: n/a
  description: Return value of `_partition_needs_state_sample`.
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

# _partition_needs_state_sample

## Purpose
Tests whether any effector in a solver partition uses the typed wrench interface and therefore needs a `StateSample` built before evaluation.

## Design & Implementation
An `any` over the effector tuple requiring both `_effector_in_partition` and `_wrench_method_available`. Declared `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `partition` | Symbol | n/a | yes | Positional argument `partition`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_partition_needs_state_sample`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1026-1026`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:128-128`

**Downstream**

- `callees` → [[simulation.effector_sampling__wrench_method_available|_wrench_method_available]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:39-39`
<!-- vulcan:connections:end -->

## Limitations
`_wrench_method_available` performs a `hasmethod` query, which is not free on every RHS call.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 35.
