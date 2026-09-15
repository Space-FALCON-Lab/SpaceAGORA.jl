---
id: gnc.mesh_distance__rpo_station_point
label: _rpo_station_point
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/mesh_distance.jl
  symbol: _rpo_station_point
  lines:
  - 2
  - 2
inputs:
- id: points
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `points`.
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
  type: SVector
  units: n/a
  description: Return value of `_rpo_station_point`. Returns `SVector{3, Float64}(points[1,
    idx], points[2, idx], points[3, idx])`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# _rpo_station_point

## Purpose
Extracts a single station point from the column-major `3 x N` point-cloud matrix as a stack-allocated `SVector{3, Float64}`, so the KD-tree searchers can do vector arithmetic without heap allocation.

## Design & Implementation
Marked `@inline` and written as `SVector{3, Float64}(points[1, idx], points[2, idx], points[3, idx])`, reading the three components of column `idx` explicitly rather than slicing with `points[:, idx]`, which would allocate a `Vector`. Callers are `_nearest_station_point_kdtree`, which uses the result in `sum(abs2, p - q)`, and `nearest_station_point`, which returns it to the caller as the winning point.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Matrix{Float64} | n/a | yes | Positional argument `points`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `_rpo_station_point`. Returns `SVector{3, Float64}(points[1, idx], points[2, idx], points[3, idx])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.mesh_distance_nearest_station_point|nearest_station_point]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:89-89`
- [[gncz.mesh_distance__nearest_station_point_kdtree|_nearest_station_point_kdtree]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:25-25`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Indexing is unchecked in the sense that no bounds validation is done here beyond Julia's default array bounds check, so an out-of-range `idx` throws a `BoundsError` with no station-level context. It assumes the matrix has exactly three rows; a matrix with more rows silently reads only the first three, and the invariant is enforced only once, at `RPOStationGeometry` construction.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl` line 2.
