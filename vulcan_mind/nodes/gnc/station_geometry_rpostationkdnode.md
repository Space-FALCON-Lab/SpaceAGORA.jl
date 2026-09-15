---
id: gnc.station_geometry_rpostationkdnode
label: RPOStationKDNode
kind: struct
source:
  file: src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl
  symbol: RPOStationKDNode
  lines:
  - 2
  - 2
inputs:
- id: idx
  type: Int
  units: n/a
  required: true
  description: Field `idx`.
- id: axis
  type: Int
  units: n/a
  required: true
  description: Field `axis`.
- id: left
  type: Union{Nothing, RPOStationKDNode}
  units: n/a
  required: true
  description: Field `left`.
- id: right
  type: Union{Nothing, RPOStationKDNode}
  units: n/a
  required: true
  description: Field `right`.
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
  description: Constructed `RPOStationKDNode`.
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

# RPOStationKDNode

## Purpose
Node type of the balanced KD-tree that indexes a station's body-frame point cloud, enabling logarithmic nearest-point queries used by RPO keep-out and proximity checks instead of a linear scan over every sample.

## Design & Implementation
A `mutable struct` with four fields: `idx::Int`, the column index of this node's point in the owning `3 x N` `points_body` matrix; `axis::Int`, the splitting dimension in `1:3`; and `left` / `right`, each `Union{Nothing, RPOStationKDNode}` so that leaves terminate with `nothing` rather than a sentinel node. Storing only an index keeps the tree free of duplicated coordinate data, and the union-typed children let the recursive searchers in `mesh_distance.jl` dispatch on `::Nothing` for the base case.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `idx` | Int | n/a | yes | Field `idx`. |
| in | `axis` | Int | n/a | yes | Field `axis`. |
| in | `left` | Union{Nothing, RPOStationKDNode} | n/a | yes | Field `left`. |
| in | `right` | Union{Nothing, RPOStationKDNode} | n/a | yes | Field `right`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOStationKDNode | n/a | — | Constructed `RPOStationKDNode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.station_geometry__rpo_build_station_kdtree|_rpo_build_station_kdtree]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl:25-25`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct is mutable but exposes no rebalancing or insertion operation, so the tree is effectively immutable after `_rpo_build_station_kdtree` returns; adding points requires a full rebuild. The `Union{Nothing, RPOStationKDNode}` children force a pointer-chasing layout with branch-unpredictable recursion, and nothing validates that `axis` is in `1:3` or that `idx` is a legal column of the matrix it will later be used to index.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl` line 2.
