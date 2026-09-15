---
id: simulation.dynamics_rhs__rhs_flat_item_estimated_cost_ns
label: _rhs_flat_item_estimated_cost_ns
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_flat_item_estimated_cost_ns
  lines:
  - 451
  - 451
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: item
  type: Int
  units: n/a
  required: true
  description: Positional argument `item`.
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
  type: Float64
  units: n/a
  description: Return value of `_rhs_flat_item_estimated_cost_ns`.
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

# _rhs_flat_item_estimated_cost_ns

## Purpose
Estimated cost of one flat work item, derived from its effector's estimate since the cost model is per effector rather than per satellite.

## Design & Implementation
Extracts the effector index from the item and calls `_rhs_effector_estimated_cost_ns`. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `item` | Int | n/a | yes | Positional argument `item`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_rhs_flat_item_estimated_cost_ns`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__prepare_rhs_flat_work_packets_bang|_prepare_rhs_flat_work_packets!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:516-516`
- [[simulation.dynamics_rhs__rhs_flat_packet_work_stats|_rhs_flat_packet_work_stats]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:556-556`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__rhs_effector_estimated_cost_ns|_rhs_effector_estimated_cost_ns]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:453-453`
- `callees` → [[simulation.dynamics_rhs__rhs_flat_item_eff_idx|_rhs_flat_item_eff_idx]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:452-452`
<!-- vulcan:connections:end -->

## Limitations
Ignores per-satellite variation such as one satellite being in the atmosphere and another not.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 451.
