---
id: simulation.dynamics_rhs__rhs_flat_item_eff_idx
label: _rhs_flat_item_eff_idx
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_flat_item_eff_idx
  lines:
  - 447
  - 447
inputs:
- id: item
  type: Int
  units: n/a
  required: true
  description: Positional argument `item`.
- id: n_effectors
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_effectors`.
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
  description: Return value of `_rhs_flat_item_eff_idx`.
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

# _rhs_flat_item_eff_idx

## Purpose
Recovers the effector index from a flat work item within the cost-model code, kept as a named alias so that code reads in its own vocabulary rather than the constellation-node one.

## Design & Implementation
Forwards `item` and `n_effectors` to `_constellation_node_eff_idx`, which returns `mod(item - 1, n_effectors) + 1`. Declared `@inline` with an `::Int` return, so the alias has no runtime cost.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `item` | Int | n/a | yes | Positional argument `item`. |
| in | `n_effectors` | Int | n/a | yes | Positional argument `n_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_flat_item_eff_idx`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__rhs_flat_item_estimated_cost_ns|_rhs_flat_item_estimated_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:452-452`
- [[simulation.dynamics_rhs__update_rhs_flat_packet_cost_model_bang|_update_rhs_flat_packet_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:636-636`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__constellation_node_eff_idx|_constellation_node_eff_idx]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:448-448`
<!-- vulcan:connections:end -->

## Limitations
Two names for one decode make the satellite-major layout contract slightly harder to trace when reading the packet cost-model code.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 447.
