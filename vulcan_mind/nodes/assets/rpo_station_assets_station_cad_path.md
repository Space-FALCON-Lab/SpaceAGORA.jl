---
id: assets.rpo_station_assets_station_cad_path
label: station_cad_path
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: station_cad_path
  lines:
  - 20
  - 20
inputs:
- id: kind
  type: Symbol
  units: n/a
  required: false
  description: Positional argument `kind` (default `:gateway`).
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
  description: Return value of `station_cad_path`. Returns `joinpath(_STATION_GEOMETRY_ROOT,
    "gateway", "Gateway_Core.stl")`.
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

# station_cad_path

## Purpose
Resolves the on-disk STL file for a CAD-backed station kind.

## Design & Implementation
Defaults `kind` to `:gateway`, accepts `:gateway_core` as an alias, and returns `gateway/Gateway_Core.stl` under the geometry root for either. Everything else raises `ArgumentError`. Keeping the alias here rather than in callers means the two names cannot drift to different files.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kind` | Symbol | n/a | no | Positional argument `kind` (default `:gateway`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `station_cad_path`. Returns `joinpath(_STATION_GEOMETRY_ROOT, "gateway", "Gateway_Core.stl")`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[assets.rpo_station_assets_load_rpo_station_cad_pointcloud|load_rpo_station_cad_pointcloud]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:161-161`
- [[assets.rpo_station_assets_load_rpo_station_cad_triangles|load_rpo_station_cad_triangles]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:152-152`
- [[assets.rpo_station_assets_station_geometry_path|station_geometry_path]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:15-15`
- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Existence is not checked here — the loaders check `isfile` themselves — so a caller that uses the path for anything other than the provided loaders must test for the file.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 20.
