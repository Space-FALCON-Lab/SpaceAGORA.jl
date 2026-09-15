---
id: assets.rpo_station_assets__load_stl_triangles
label: _load_stl_triangles
kind: function
source:
  file: src/assets/rpo_station_assets.jl
  symbol: _load_stl_triangles
  lines:
  - 54
  - 54
inputs:
- id: path
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `path`.
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
  description: Return value of `_load_stl_triangles`. Returns `open(path, "r") do
    io` or `triangles`.
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

# _load_stl_triangles

## Purpose
Parses an STL file in either encoding into a three-by-3N matrix of vertex columns, three consecutive columns per triangle.

## Design & Implementation
If `_stl_is_binary` says binary, it skips the 80-byte header, reads the `UInt32` count, and for each triangle reads a `Float32` normal it discards, three `Float32` vertices it widens into the output columns, and a `UInt16` attribute it discards. Otherwise it scans lines, strips them, and for each one beginning with `vertex` parses the next three whitespace-separated tokens as `Float64`, collecting vertices in order and copying them into the matrix at the end. Both paths produce the same column layout so downstream code is encoding-agnostic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | AbstractString | n/a | yes | Positional argument `path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_load_stl_triangles`. Returns `open(path, "r") do io` or `triangles`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[assets.rpo_station_assets__stl_pointcloud|_stl_pointcloud]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:118-118`
- [[assets.rpo_station_assets_load_rpo_station_cad_triangles|load_rpo_station_cad_triangles]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:154-154`
- [[module.assets|RPOStationAssets]] · `api` → `module_api` · call · `src/assets/rpo_station_assets.jl`

**Downstream**

- `callees` → [[assets.rpo_station_assets__stl_is_binary|_stl_is_binary]] · `callers` · call · `src/assets/rpo_station_assets.jl:55-55`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/assets/rpo_station_assets.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
The ASCII path trusts that vertices appear in groups of three inside `facet` blocks and does not verify the count is divisible by three; the binary path allocates three small `Float32` vectors per triangle rather than reusing buffers, which is slow for large meshes. Stored normals are discarded, so winding order is the only orientation information retained.

## Provenance
Mapped from `src/assets/rpo_station_assets.jl` line 54.
