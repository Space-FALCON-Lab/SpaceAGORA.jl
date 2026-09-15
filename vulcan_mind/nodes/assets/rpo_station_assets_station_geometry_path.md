---
id: assets.rpo_station_assets_station_geometry_path
label: station_geometry_path
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: station_geometry_path
  lines:
  - 11
  - 11
inputs:
- id: kind
  type: Symbol
  units: n/a
  required: false
  description: Positional argument `kind` (default `:demo`).
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
  description: Return value of `station_geometry_path`. Returns `joinpath(_STATION_GEOMETRY_ROOT,
    "demo", "station_pointcloud.csv")` or `station_cad_path(kind)`.
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

# station_geometry_path

## Purpose
Maps a station kind symbol onto the file that holds its geometry, whichever representation that station ships in.

## Design & Implementation
Defaults `kind` to `:demo` and returns the CSV point cloud path under `demo/`. For `:gateway` or `:gateway_core` it delegates to `station_cad_path`, so the same symbol resolves regardless of whether the caller asked through the geometry or the CAD accessor. Any other symbol raises `ArgumentError` naming the valid kinds and pointing at the artifact-loader extension route.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kind` | Symbol | n/a | no | Positional argument `kind` (default `:demo`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `station_geometry_path`. Returns `joinpath(_STATION_GEOMETRY_ROOT, "demo", "station_pointcloud.csv")` or `station_cad_path(kind)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[assets.load_rpo_station_pointcloud|load_rpo_station_pointcloud]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:29-29`

**Downstream**

- `callees` → [[assets.rpo_station_assets_station_cad_path|station_cad_path]] · `callers` · call · `src/assets/rpo_station_assets.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
The return type differs by kind — a CSV for demo, an STL for gateway — and nothing in the return value indicates which, so a caller must already know the format to parse the file it is handed.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 11.
