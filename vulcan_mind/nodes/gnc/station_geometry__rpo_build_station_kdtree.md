---
id: gnc.station_geometry__rpo_build_station_kdtree
label: _rpo_build_station_kdtree
kind: function
source:
  file: src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl
  symbol: _rpo_build_station_kdtree
  lines:
  - 18
  - 18
inputs:
- id: points
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `points`.
- id: indices
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `indices`.
- id: depth
  type: Int
  units: n/a
  required: false
  description: Positional argument `depth` (default `0`).
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
  type: RPOStationKDNode
  units: n/a
  description: Return value of `_rpo_build_station_kdtree`. Returns `RPOStationKDNode(`.
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

# _rpo_build_station_kdtree

## Purpose
Builds the KD-tree over a station point cloud by recursively median-splitting the supplied column indices, producing the `kd_root` stored in `RPOStationGeometry` for nearest-point queries.

## Theory & Math
Splitting axis at depth $d$ is $a = (d \bmod 3) + 1$, and the pivot is the element at rank $m = \lceil n/2 \rceil$ of the indices sorted by coordinate $x_a$, where $n$ is the number of indices in the current subtree. This yields tree height $O(\log_2 N)$ for $N$ points.

## Design & Implementation
Takes `points::Matrix{Float64}` (3 x N), a `Vector{Int}` of column indices, and a `depth::Int` defaulting to 0. Empty `indices` returns `nothing`, terminating a branch. The split axis cycles through the three spatial dimensions as `mod(depth, 3) + 1`. It calls `sort!(indices; by = idx -> points[axis, idx])`, mutating the caller's vector in place, then picks `mid = cld(length(indices), 2)` as the pivot and recurses on the slices `indices[1:mid-1]` and `indices[mid+1:end]`, each at `depth + 1`. The top-level caller in `RPOStationGeometry` passes `collect(1:size(points, 2))`, a fresh vector, so the in-place sorting is safe there.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Matrix{Float64} | n/a | yes | Positional argument `points`. |
| in | `indices` | Vector{Int} | n/a | yes | Positional argument `indices`. |
| in | `depth` | Int | n/a | no | Positional argument `depth` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOStationKDNode | n/a | — | Return value of `_rpo_build_station_kdtree`. Returns `RPOStationKDNode(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl:38-38`
- `callees` → [[gnc.station_geometry_rpostationkdnode|RPOStationKDNode]] · `callers` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl:25-25`
- `callees` → [[gncz.station_geometry_rpostationgeometry|RPOStationGeometry]] · `callers` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl:34-34`
<!-- vulcan:connections:end -->

## Limitations
`sort!` mutates the index vector the caller passed, so calling this directly with a vector you still need will scramble its order. Each recursion level allocates two new index slices, giving O(N log N) temporary allocation on top of the sort, and the recursion depth is proportional to log N but unbounded by any explicit guard. The cyclic `mod(depth, 3)` axis choice ignores per-axis variance, so a point cloud that is nearly planar or strongly elongated in one axis produces poorly balanced splitting planes and slower queries than a variance-selected axis would.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl` line 18.
