---
id: gncz.mesh_distance__nearest_station_point_kdtree
label: _nearest_station_point_kdtree
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/mesh_distance.jl
  symbol: _nearest_station_point_kdtree
  lines:
  - 7
  - 41
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks namespace providing the station geometry, its KD-tree
    root, and the point-cloud storage matrix.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: nearest
  type: Tuple{Int, Float64}
  units: index, m^2
  description: Index of the closest station point-cloud sample and the squared distance
    from the query point to it.
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

# _nearest_station_point_kdtree

## Purpose
`_nearest_station_point_kdtree` performs the branch-and-bound descent that turns the station point cloud into a fast nearest-neighbour query. Every clearance, surface-normal, and standoff computation in the proximity-operations navigation stack bottoms out here, so its cost dominates path scoring when a planner evaluates thousands of candidate samples.

## Theory & Math
The search maintains the best index and best squared distance found so far. At each node it evaluates the squared distance to the stored point, descends first into the half-space containing the query along the node splitting axis, and only visits the far half-space when the squared axis offset $(p_a - q_a)^2$ is smaller than the current best squared distance. That test is exactly the condition under which the far side can still contain a closer point, so the pruned search is exact rather than approximate.

## Model & Assumptions
Two methods are defined, one dispatching on a missing child to terminate the recursion and one on a real node, which removes any branch on emptiness from the hot path. Queries and stored points are three-dimensional and the tree is assumed to have been built over the same point matrix that is passed in. Distances stay squared throughout so no square root is taken until the caller needs a metric distance.

## Design & Implementation
A sibling pair of methods computes the squared distance alone without tracking the index, avoiding the static-vector fetch when only a margin is needed. The two public wrappers convert the query to a static three-vector, reject an empty tree with an argument error, and take the square root once at the end.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks namespace providing the station geometry, its KD-tree root, and the point-cloud storage matrix. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `nearest` | Tuple{Int, Float64} | index, m^2 | — | Index of the closest station point-cloud sample and the squared distance from the query point to it. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.mesh_distance_nearest_station_point|nearest_station_point]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:88-88`

**Downstream**

- `callees` → [[gnc.mesh_distance__rpo_station_point|_rpo_station_point]] · `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:25-25`
<!-- vulcan:connections:end -->

## Limitations
The recursion is not tail-call eliminated, so a badly unbalanced tree can deepen the stack. The tree is built once from a fixed point matrix and cannot absorb a moving or deforming station without a rebuild.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:1-91`.
