---
id: gnc.mesh_distance_nearest_station_point
label: nearest_station_point
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/mesh_distance.jl
  symbol: nearest_station_point
  lines:
  - 85
  - 85
inputs:
- id: p_body
  type: Any
  units: n/a
  required: true
  description: Positional argument `p_body`.
- id: station
  type: RPOStationGeometry
  units: n/a
  required: true
  description: Positional argument `station`.
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
  description: Return value of `nearest_station_point`. Returns `q, sqrt(best_d2),
    best_idx`.
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

# nearest_station_point

## Purpose
Returns the full nearest-neighbour result for a body-frame query against a station point cloud: the closest point itself, its true distance, and its column index, for callers that need the contact geometry rather than just a range.

## Design & Implementation
Converts `p_body` to `SVector{3, Float64}` and throws `ArgumentError("RPO station KD-tree is empty.")` if `station.kd_root === nothing`. It calls `_nearest_station_point_kdtree` seeded with `best_idx = 1` and `best_d2 = Inf`, re-fetches the winning coordinates via `_rpo_station_point`, and returns the tuple `(q, sqrt(best_d2), best_idx)`. The `sqrt` is applied once, outside the recursion, so the search itself stays in squared-distance space.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p_body` | Any | n/a | yes | Positional argument `p_body`. |
| in | `station` | RPOStationGeometry | n/a | yes | Positional argument `station`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `nearest_station_point`. Returns `q, sqrt(best_d2), best_idx`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.clearance_rpo_clearance_to_station|rpo_clearance_to_station]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:3-3`
- [[gncz.surface_frames_rpo_surface_normal_from_pointcloud|rpo_surface_normal_from_pointcloud]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/surface_frames.jl:3-3`

**Downstream**

- `callees` → [[gnc.mesh_distance__rpo_station_point|_rpo_station_point]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:89-89`
- `callees` → [[gncz.mesh_distance__nearest_station_point_kdtree|_nearest_station_point_kdtree]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
The incumbent index is seeded to `1` rather than a sentinel, so if the tree traversal never improves on `Inf` the function silently reports column 1 as the answer; this is reachable when the query contains NaN, since all distance comparisons then evaluate false. The returned index refers into `station.points_body` and becomes meaningless if the geometry is rebuilt. As with the distance-only variant, `keepout_radius_m` is ignored and the answer is only as accurate as the sampling density of the point cloud.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl` line 85.
