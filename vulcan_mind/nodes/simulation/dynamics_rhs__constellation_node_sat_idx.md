---
id: simulation.dynamics_rhs__constellation_node_sat_idx
label: _constellation_node_sat_idx
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _constellation_node_sat_idx
  lines:
  - 388
  - 388
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
  description: Return value of `_constellation_node_sat_idx`.
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

# _constellation_node_sat_idx

## Purpose
Recovers the satellite index from a satellite-major flat work-item number, the inverse of `_constellation_node_work_item` for the satellite component.

## Design & Implementation
Returns `fld(item - 1, n_effectors) + 1`. Declared `@inline` with an `::Int` return. Used by the flat queue worker to find which satellite's partial buffer an item's result belongs in.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `item` | Int | n/a | yes | Positional argument `item`. |
| in | `n_effectors` | Int | n/a | yes | Positional argument `n_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_constellation_node_sat_idx`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1104-1104`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Shares the layout coupling of its siblings.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 388.
