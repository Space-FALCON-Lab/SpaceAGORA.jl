---
id: simulation.dynamics_rhs__flat_partition_selected
label: _flat_partition_selected
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _flat_partition_selected
  lines:
  - 325
  - 325
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
- id: partition
  type: Union{Nothing, Symbol}
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
  description: Return value of `_flat_partition_selected`.
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

# _flat_partition_selected

## Purpose
Tests whether an effector should be included in the current flat-queue RHS call, honouring the split solver's implicit-or-explicit partition when one is given and selecting every effector when the call is the monolithic slow RHS.

## Design & Implementation
Returns true immediately for a `nothing` partition, otherwise delegates to `_effector_in_partition(effector, partition)`, which compares the effector's validated `solver_partition` against the requested symbol. Declared `@inline` with a `::Bool` return so the work-item builder's inner loop pays nothing for the check.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `partition` | Union{Nothing, Symbol} | n/a | yes | Positional argument `partition`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_flat_partition_selected`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__prepare_rhs_flat_work_items_bang|_prepare_rhs_flat_work_items!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:414-414`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the `nothing` case short-circuits before validation, an effector with a mis-declared partition is only caught when the split solver is actually in use, not in the ordinary monolithic path.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 325.
