---
id: gnc.mesh_distance_nearest_station_distance_sq
label: nearest_station_distance_sq
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/mesh_distance.jl
  symbol: nearest_station_distance_sq
  lines:
  - 78
  - 78
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
  description: Return value of `nearest_station_distance_sq`. Returns `_nearest_station_distance_sq_kdtree(station.kd_root,
    p, station.points_body, Inf`.
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

# nearest_station_distance_sq

## Purpose
Public entry point returning the squared distance in metres-squared from a body-frame query point to the closest point of an `RPOStationGeometry` point cloud, the quantity RPO collision and keep-out predicates compare against a squared threshold.

## Design & Implementation
Converts `p_body` to `SVector{3, Float64}`, so any 3-element iterable is accepted. It then asserts the tree exists, throwing `ArgumentError("RPO station KD-tree is empty.")` when `station.kd_root === nothing`, and delegates to `_nearest_station_distance_sq_kdtree(station.kd_root, p, station.points_body, Inf)` with `Inf` as the initial incumbent. Returning the squared distance rather than the distance deliberately avoids a `sqrt` in the inner loop of callers that only threshold the result.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p_body` | Any | n/a | yes | Positional argument `p_body`. |
| in | `station` | RPOStationGeometry | n/a | yes | Positional argument `station`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `nearest_station_distance_sq`. Returns `_nearest_station_distance_sq_kdtree(station.kd_root, p, station.points_body, Inf`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:10-10`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl`

**Downstream**

- `callees` → [[gnc.mesh_distance__nearest_station_distance_sq_kdtree|_nearest_station_distance_sq_kdtree]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:81-81`
<!-- vulcan:connections:end -->

## Limitations
The station's `keepout_radius_m` field is not applied here; callers must add or compare it themselves, which is easy to forget. Converting `p_body` with the `SVector{3, Float64}` constructor throws if the input does not have exactly three elements, and a query containing NaN propagates as `Inf` rather than an error because every comparison in the search kernel fails. Only the point samples are considered, so concave mesh regions between samples are invisible and the reported distance is an upper bound on the true surface distance.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl` line 78.
