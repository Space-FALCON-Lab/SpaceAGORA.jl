---
id: gncz.station_geometry_rpostationgeometry
label: RPOStationGeometry
kind: struct
source:
  file: src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl
  symbol: RPOStationGeometry
  lines:
  - 10
  - 15
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks namespace where the station point cloud type and its
    KD-tree node type are declared.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: station_geometry
  type: RPOStationGeometry
  units: m
  description: Validated station point cloud, its balanced KD-tree root, the keepout
    radius, and the station name.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# RPOStationGeometry

## Purpose
`RPOStationGeometry` is the target-side collision model for proximity operations. It owns the body-frame point cloud that samples the station surface, the spatial index that makes nearest-point queries affordable, and the keepout radius that inflates the structure into a volume no chaser may enter.

## Theory & Math
The index is a median-split KD-tree. At depth $d$ the splitting axis is $(d \bmod 3) + 1$, the surviving indices are sorted by that coordinate, and the median index becomes the node while the two halves recurse with depth increased by one. Median splitting gives a tree of depth logarithmic in the point count, which is what bounds the branch-and-bound nearest-point search.

## Model & Assumptions
Points are stored as a three-by-N matrix of double precision coordinates in the station body frame. The keyword constructor validates that the matrix has exactly three rows, that it contains at least one point, and that the keepout radius is non-negative, raising an argument error otherwise so a malformed cloud cannot reach the distance loops. The KD-tree root may be nothing only for a cloud that failed to build, and the distance queries reject that case explicitly.

## Design & Implementation
The mutable KD-tree node stores a point index, the splitting axis, and two children typed as either a node or nothing, so the recursion terminates by dispatch. Building sorts index vectors rather than moving points, leaving the coordinate matrix contiguous for cache-friendly distance evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks namespace where the station point cloud type and its KD-tree node type are declared. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `station_geometry` | RPOStationGeometry | m | — | Validated station point cloud, its balanced KD-tree root, the keepout radius, and the station name. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning_rpo_geometry_with_replanning_spheres|rpo_geometry_with_replanning_spheres]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:163-163`
- [[gnc.station_geometry__rpo_build_station_kdtree|_rpo_build_station_kdtree]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl:34-34`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Sorting at every level makes construction more expensive than a linear-time median selection, and the tree is immutable in practice because no insertion or rebalancing routine exists. The point cloud is static, so a rotating or articulated station cannot be represented without rebuilding the whole geometry.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl:1-42`.
