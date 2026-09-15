---
id: simulation.dynamics_rhs__prepare_rhs_flat_work_items_bang
label: _prepare_rhs_flat_work_items!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _prepare_rhs_flat_work_items!
  lines:
  - 396
  - 396
inputs:
- id: work_items
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `work_items`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: Int
  units: n/a
  description: Return value of `_prepare_rhs_flat_work_items!`; mutates `work_items`
    in place.
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

# _prepare_rhs_flat_work_items!

## Purpose
Builds the list of satellite-effector work items for the flat queue, excluding inactive satellites and effectors handled by batch kernels.

## Design & Implementation
Resizes the item vector, loops active satellites and selected effectors, skips batchable and harmonics-prepass effectors when no partition is given, and encodes items. Returns the count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `work_items` | Vector{Int} | n/a | yes | Positional argument `work_items`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `partition` | Union{Nothing, Symbol} | n/a | yes | Positional argument `partition`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_prepare_rhs_flat_work_items!`; mutates `work_items` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__build_constellation_execution_plan_bang|_build_constellation_execution_plan!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:465-465`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__constellation_node_eff_idx|_constellation_node_eff_idx]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:440-440`
- `callees` → [[simulation.dynamics_rhs__constellation_node_work_item|_constellation_node_work_item]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:418-418`
- `callees` → [[simulation.dynamics_rhs__count_flat_queue_only_effectors|_count_flat_queue_only_effectors]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:429-429`
- `callees` → [[simulation.dynamics_rhs__flat_partition_selected|_flat_partition_selected]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:414-414`
- `callees` → [[simulation.dynamics_rhs__rhs_effector_estimated_cost_ns|_rhs_effector_estimated_cost_ns]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:435-435`
<!-- vulcan:connections:end -->

## Limitations
Item order is satellite-major, which affects packet grouping.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 396.
