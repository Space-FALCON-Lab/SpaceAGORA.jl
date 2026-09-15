---
id: assets.rpo_station_assets_load_rpo_station_cad_pointcloud
label: load_rpo_station_cad_pointcloud
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: load_rpo_station_cad_pointcloud
  lines:
  - 160
  - 160
inputs:
- id: kind
  type: Symbol
  units: n/a
  required: false
  description: Positional argument `kind` (default `:gateway`).
- id: n_points
  type: Integer
  units: n/a
  required: false
  description: Keyword argument `n_points` (default `10000`).
- id: rng
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rng` (default `MersenneTwister(42)`).
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
  description: Return value of `load_rpo_station_cad_pointcloud`. Returns `_stl_pointcloud(path;
    n_points=n_points, rng=rng, center=center, scale=scale)`.
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

# load_rpo_station_cad_pointcloud

## Purpose
Public loader producing a reproducible surface point cloud of a CAD-backed station for clearance and visualisation.

## Design & Implementation
Resolves and existence-checks the STL path, then delegates to `_stl_pointcloud` forwarding `n_points`, `rng`, `center` and `scale`. The `rng` default is `MersenneTwister(42)`, a fixed seed, so two runs of the same scenario see the identical cloud and clearance statistics are comparable across runs; a caller wanting fresh samples passes its own generator.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kind` | Symbol | n/a | no | Positional argument `kind` (default `:gateway`). |
| in | `n_points` | Integer | n/a | no | Keyword argument `n_points` (default `10000`). |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `MersenneTwister(42)`). |
| in | `center` | Bool | n/a | no | Keyword argument `center` (default `true`). |
| in | `scale` | Real | n/a | no | Keyword argument `scale` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `load_rpo_station_cad_pointcloud`. Returns `_stl_pointcloud(path; n_points=n_points, rng=rng, center=center, scale=scale)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- `callees` → [[assets.rpo_station_assets__stl_pointcloud|_stl_pointcloud]] · `callers` · call · `src/assets/rpo_station_assets.jl:163-163`
- `callees` → [[assets.rpo_station_assets_station_cad_path|station_cad_path]] · `callers` · call · `src/assets/rpo_station_assets.jl:161-161`
<!-- vulcan:connections:end -->

## Limitations
The default of 10,000 points is a fixed budget unrelated to the mesh's size, so a large station is sampled sparsely and a small one densely; the fixed seed means independent Monte Carlo runs share the same cloud unless the caller remembers to override it.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 160.
