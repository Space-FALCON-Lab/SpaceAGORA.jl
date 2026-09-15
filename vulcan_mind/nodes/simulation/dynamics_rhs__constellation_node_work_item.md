---
id: simulation.dynamics_rhs__constellation_node_work_item
label: _constellation_node_work_item
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _constellation_node_work_item
  lines:
  - 384
  - 384
inputs:
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: eff_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `eff_idx`.
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
  description: Return value of `_constellation_node_work_item`.
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

# _constellation_node_work_item

## Purpose
Encodes a satellite-effector pair as one integer work item, satellite-major, so the flat queue can carry its work list as a plain `Vector{Int}`.

## Design & Implementation
Returns `(sat_idx - 1) * n_effectors + eff_idx`. Satellite-major ordering keeps all of one satellite's effectors adjacent, which the packet grouper exploits so a packet tends to touch few satellites' state. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `eff_idx` | Int | n/a | yes | Positional argument `eff_idx`. |
| in | `n_effectors` | Int | n/a | yes | Positional argument `n_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_constellation_node_work_item`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__prepare_rhs_flat_work_items_bang|_prepare_rhs_flat_work_items!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:418-418`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The maximum item is `num_sats * n_effectors`, comfortably within `Int` but the layout is implicit and shared with the two decoders.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 384.
