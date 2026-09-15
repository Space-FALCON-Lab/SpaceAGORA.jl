---
id: gnc.mesh_distance__nearest_station_distance_sq_kdtree
label: _nearest_station_distance_sq_kdtree
kind: function
source:
  file: src/gnc/navigation/rpo_nav/distances/mesh_distance.jl
  symbol: _nearest_station_distance_sq_kdtree
  lines:
  - 44
  - 44
inputs:
- id: node
  type: Nothing
  units: n/a
  required: true
  description: Positional argument `node`.
- id: p
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `p`.
- id: points
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Positional argument `points`.
- id: best_d2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `best_d2`.
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
  description: Return value of `_nearest_station_distance_sq_kdtree`. Returns `best_d2`.
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

# _nearest_station_distance_sq_kdtree

## Purpose
Recursive KD-tree search kernel that returns the smallest squared Euclidean distance from a query point to any point in the station cloud, used by proximity and keep-out checks that never need the identity of the winning point.

## Theory & Math
Pruning is justified by the fact that every point in the far subtree lies on the opposite side of the splitting plane $x_a = q_a$, so its squared distance is at least $(p_a - q_a)^2$. The far branch can therefore be skipped whenever $(p_a - q_a)^2 \ge d^2_{\text{best}}$, where $p$ is the query, $q$ the node point, and $a$ the node's split axis.

## Design & Implementation
Two methods implement the recursion. The `::Nothing` method is the base case and simply returns `best_d2` unchanged. The `::RPOStationKDNode` method computes the squared distance to the node's own point with explicit `dx`, `dy`, `dz` scalars rather than `SVector` subtraction, updates `best_d2` if closer, then computes `diff = p[axis] - points[axis, node.idx]` and descends the near branch first (`left` when `diff <= 0.0`, else `right`). The far branch is visited only when `diff * diff < best_d2`, the standard hyperplane-pruning test. `best_d2` threads through the recursion as an accumulator rather than being stored anywhere.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `node` | Nothing | n/a | yes | Positional argument `node`. |
| in | `p` | SVector{3, Float64} | n/a | yes | Positional argument `p`. |
| in | `points` | Matrix{Float64} | n/a | yes | Positional argument `points`. |
| in | `best_d2` | Float64 | n/a | yes | Positional argument `best_d2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_nearest_station_distance_sq_kdtree`. Returns `best_d2`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.mesh_distance_nearest_station_distance_sq|nearest_station_distance_sq]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl:81-81`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Recursion is unbounded by any depth guard, so a pathologically deep tree could exhaust the stack, and Julia does not perform tail-call elimination on the non-tail calls here. The scalar-per-component distance loop is faster than the `SVector` version used by `_nearest_station_point_kdtree` but duplicates that logic, so the two can drift apart. A NaN in the query point makes every comparison false, silently returning the incoming `best_d2` (`Inf` from the entry point) instead of raising.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl` line 44.
