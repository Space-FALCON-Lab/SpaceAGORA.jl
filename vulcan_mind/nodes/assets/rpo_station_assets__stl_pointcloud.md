---
id: assets.rpo_station_assets__stl_pointcloud
label: _stl_pointcloud
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: _stl_pointcloud
  lines:
  - 116
  - 116
inputs:
- id: path
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `path`.
- id: n_points
  type: Integer
  units: n/a
  required: false
  description: Keyword argument `n_points` (default `10000`).
- id: rng
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rng` (default `Random.default_rng()`).
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
  description: Return value of `_stl_pointcloud`. Returns `points`.
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

# _stl_pointcloud

## Purpose
Converts a triangle mesh into a fixed-size point cloud whose density is proportional to surface area, the representation the RPO clearance checks consume.

## Design & Implementation
Rejects a non-positive `n_points`, loads triangles, applies `scale` when it differs from one and centres when requested. It computes each triangle's area as half the norm of the cross product of two edges, rejects a mesh with no triangles or zero total area, and builds a normalised cumulative distribution. For each output point it draws a uniform, finds the triangle by `searchsortedfirst` on the CDF, clamps the index, and samples a point on that triangle. The `rng` argument makes the cloud reproducible.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | AbstractString | n/a | yes | Positional argument `path`. |
| in | `n_points` | Integer | n/a | no | Keyword argument `n_points` (default `10000`). |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `Random.default_rng()`). |
| in | `center` | Bool | n/a | no | Keyword argument `center` (default `true`). |
| in | `scale` | Real | n/a | no | Keyword argument `scale` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_stl_pointcloud`. Returns `points`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[assets.rpo_station_assets_load_rpo_station_cad_pointcloud|load_rpo_station_cad_pointcloud]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:163-163`
- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/assets/rpo_station_assets.jl:119-119`
- `callees` → [[assets.rpo_station_assets__center_triangles_bang|_center_triangles!]] · `callers` · call · `src/assets/rpo_station_assets.jl:120-120`
- `callees` → [[assets.rpo_station_assets__load_stl_triangles|_load_stl_triangles]] · `callers` · call · `src/assets/rpo_station_assets.jl:118-118`
- `callees` → [[assets.rpo_station_assets__sample_triangle_point|_sample_triangle_point]] · `callers` · call · `src/assets/rpo_station_assets.jl:141-141`
<!-- vulcan:connections:end -->

## Limitations
Degenerate triangles with zero area are kept in the CDF and simply never selected, which is correct but means the triangle count reported elsewhere overstates usable geometry; per-point column copies and broadcasts make this allocation-heavy for large `n_points`.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 116.
