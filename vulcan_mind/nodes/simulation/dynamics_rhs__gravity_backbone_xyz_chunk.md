---
id: simulation.dynamics_rhs__gravity_backbone_xyz_chunk
label: _gravity_backbone_xyz_chunk
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _gravity_backbone_xyz_chunk
  lines:
  - 1449
  - 1449
inputs:
- id: state
  type: Any
  units: n/a
  required: true
  description: Positional argument `state`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_gravity_backbone_xyz_chunk`.
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

# _gravity_backbone_xyz_chunk

## Purpose
Extracts a satellite's three-vector from whichever layout the split solver's state arrays take, since different solver configurations hand the backbone different containers.

## Design & Implementation
Handles a component tree with an `sc` field, a vector of per-satellite states exposing `pos` or `vel`, and a flat vector indexed at `3(sat_idx - 1)`, in that order. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | Any | n/a | yes | Positional argument `state`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_gravity_backbone_xyz_chunk`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__gravity_backbone_state_sample|_gravity_backbone_state_sample]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1475-1475`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Several `hasproperty` checks per call on a hot path; the flat-vector fallback assumes a three-per-satellite layout with no validation.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1449.
