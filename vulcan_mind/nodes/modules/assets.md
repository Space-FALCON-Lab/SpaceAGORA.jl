---
id: module.assets
label: RPOStationAssets
kind: module
source:
  file: src/assets/rpo_station_assets.jl
  symbol: RPOStationAssets
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: Station path helpers, point-cloud loaders, STL triangle loaders, and
    centered/scaled sampling functions exported by RPOStationAssets.
tags:
- module
charts:
- master
origin: agent
---

# RPOStationAssets

## Purpose
`RPOStationAssets` provides deterministic geometry access for rendezvous and proximity-operations examples. It resolves the repository asset paths for the demo station and Gateway CAD variants, loads text point clouds, parses binary or ASCII STL data, and exposes triangle and point-cloud representations for visualization or geometric targeting. The implementation keeps asset discovery in one module so simulation and plotting callers do not duplicate repository-relative path rules.

## Theory & Math
Triangle sampling uses barycentric area sampling. Given triangle vertices `v1`, `v2`, and `v3`, two uniform random values are transformed so the sampled point is `v1 + a*(v2-v1) + b*(v3-v1)` with the square-root correction that makes surface density uniform. Centering subtracts the mean of the loaded vertices, and scaling multiplies coordinates by the requested scalar.

## Model & Assumptions
Text point-cloud files contain three numeric coordinates per row. STL files are interpreted as triangle meshes, and binary-versus-text detection chooses the parser before triangles are materialized. Defaults are repository assets, so callers outside the checkout must supply a compatible deployment layout or call the lower-level path/loader methods with their own files.

## Design & Implementation
`station_geometry_path` and `station_cad_path` select named assets. `load_rpo_station_pointcloud` parses three-column data and validates each row. `_stl_is_binary` selects the STL reader; `_load_stl_triangles` returns a vertex matrix; `_center_triangles!` applies the requested centering. `load_rpo_station_cad_triangles` exposes the mesh, while `load_rpo_station_cad_pointcloud` calls `_stl_pointcloud` to sample triangle surfaces with a caller-provided RNG and point count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | Station path helpers, point-cloud loaders, STL triangle loaders, and centered/scaled sampling functions exported by RPOStationAssets. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[assets.rpo_station_assets__center_triangles_bang|_center_triangles!]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets__load_stl_triangles|_load_stl_triangles]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets__sample_triangle_point|_sample_triangle_point]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets__stl_is_binary|_stl_is_binary]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets__stl_pointcloud|_stl_pointcloud]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets_load_rpo_station_cad_pointcloud|load_rpo_station_cad_pointcloud]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets_load_rpo_station_cad_triangles|load_rpo_station_cad_triangles]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets_rpostationassets|RPOStationAssets]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[assets.rpo_station_assets_station_cad_path|station_cad_path]] · `module_api` · call · `src/assets/rpo_station_assets.jl`
- `api` → [[module.spaceagora|SpaceAGORA]] · `assets` · call · `src/SpaceAGORA.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
Malformed coordinate rows throw an `ArgumentError`, unsupported STL layouts can fail during parsing, and a mesh with zero-area triangles produces poor sampling statistics. Centering changes the frame origin by design; callers that need the CAD frame must pass `center=false`. The default point count of 10,000 trades memory for surface coverage and is not a physical sensor model.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl`, including the path, text-cloud, STL, centering, and sampling functions.
