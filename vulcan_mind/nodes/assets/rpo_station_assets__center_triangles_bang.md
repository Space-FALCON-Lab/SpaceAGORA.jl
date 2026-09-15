---
id: assets.rpo_station_assets__center_triangles_bang
label: _center_triangles!
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: _center_triangles!
  lines:
  - 99
  - 99
inputs:
- id: triangles
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `triangles`.
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
  description: Return value of `_center_triangles!`; mutates `triangles` in place.
    Returns `triangles`.
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

# _center_triangles!

## Purpose
Shifts a vertex matrix so the centre of its axis-aligned bounding box sits at the origin, giving station geometry a frame independent of where the CAD author placed it.

## Design & Implementation
Computes per-axis minimum and maximum across all columns with `minimum` and `maximum` along `dims=2`, forms their midpoint, and subtracts it from every column in place with a broadcast. Returns the same matrix. Bounding-box centring rather than centroid centring was chosen so the result does not depend on how finely different regions of the mesh were tessellated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `triangles` | Matrix{Float64} | n/a | yes | Positional argument `triangles`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_center_triangles!`; mutates `triangles` in place. Returns `triangles`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[assets.rpo_station_assets__stl_pointcloud|_stl_pointcloud]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:120-120`
- [[assets.rpo_station_assets_load_rpo_station_cad_triangles|load_rpo_station_cad_triangles]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:156-156`
- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The bounding-box centre is not the mass centre or the docking-port location, so a keep-out sphere defined about the origin after centring is about a geometric midpoint that may be far from any physically meaningful point on an asymmetric station.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 99.
