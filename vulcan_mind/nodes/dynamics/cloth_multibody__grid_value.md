---
id: dynamics.cloth_multibody__grid_value
label: _grid_value
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _grid_value
  lines:
  - 331
  - 331
inputs:
- id: v
  type: Any
  units: n/a
  required: true
  description: Positional argument `v`.
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  type: Any
  units: n/a
  description: Return value of `_grid_value`. Returns `v[idx]` or `v`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _grid_value

## Purpose
Lets a grid builder parameter be either one scalar for every tile or a per-tile array, indexed by linear tile index.

## Design & Implementation
Returns `v[idx]` for an `AbstractVector` or `AbstractMatrix` and `v` itself otherwise. `@inline`. The matrix branch uses linear indexing so a rows-by-cols matrix maps to the same tile order the builder uses.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v` | Any | n/a | yes | Positional argument `v`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_grid_value`. Returns `v[idx]` or `v`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_build_rectangular_compliant_grid|build_rectangular_compliant_grid]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:374-374`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Linear indexing of a matrix is column-major in Julia while the builder numbers tiles row-major, so a per-tile matrix must be transposed by the caller to line up.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 331.
