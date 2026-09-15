---
id: assets.rpo_station_assets_load_rpo_station_cad_triangles
label: load_rpo_station_cad_triangles
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: load_rpo_station_cad_triangles
  lines:
  - 151
  - 151
inputs:
- id: kind
  type: Symbol
  units: n/a
  required: false
  description: Positional argument `kind` (default `:gateway`).
- id: center
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `center` (default `true`).
- id: scale
  type: Real
  units: n/a
  required: false
  description: Keyword argument `scale` (default `1.0`).
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
  description: Return value of `load_rpo_station_cad_triangles`. Returns `triangles`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- assets
charts:
- assets
origin: agent
---

# load_rpo_station_cad_triangles

## Purpose
Public loader returning the raw triangle vertices of a CAD-backed station, for callers that need the mesh rather than a point sample.

## Design & Implementation
Resolves the STL through `station_cad_path`, raises `ArgumentError` if the file is missing, parses it with `_load_stl_triangles`, then applies `scale` and `_center_triangles!` under the same keyword defaults the point-cloud loader uses. Returning the same three-by-3N layout as the private parser means the triangle and point-cloud loaders agree on frame and units for a given `kind`, `scale` and `center`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kind` | Symbol | n/a | no | Positional argument `kind` (default `:gateway`). |
| in | `center` | Bool | n/a | no | Keyword argument `center` (default `true`). |
| in | `scale` | Real | n/a | no | Keyword argument `scale` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `load_rpo_station_cad_triangles`. Returns `triangles`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/assets/rpo_station_assets.jl:155-155`
- `callees` → [[assets.rpo_station_assets__center_triangles_bang|_center_triangles!]] · `callers` · call · `src/assets/rpo_station_assets.jl:156-156`
- `callees` → [[assets.rpo_station_assets__load_stl_triangles|_load_stl_triangles]] · `callers` · call · `src/assets/rpo_station_assets.jl:154-154`
- `callees` → [[assets.rpo_station_assets_station_cad_path|station_cad_path]] · `callers` · call · `src/assets/rpo_station_assets.jl:152-152`
<!-- vulcan:connections:end -->

## Limitations
The full mesh is re-parsed from disk on every call with no caching, so a scenario that calls this repeatedly pays the STL parse each time.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 151.
