---
id: gnc.rrt_connect_rpo_rrt_near_indices
label: rpo_rrt_near_indices
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_near_indices
  lines:
  - 49
  - 49
inputs:
- id: tree
  type: RPORRTConnectTree
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: q
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `q`.
- id: radius_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `radius_m`.
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
  description: Return value of `rpo_rrt_near_indices`. Returns `hypr_rrt_near_indices(tree,
    q, max(Float64(radius_m), 0.0))`.
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

# rpo_rrt_near_indices

## Purpose
Returns the indices of all tree nodes within `radius_m` of the query `q`, supplying the candidate parent and rewiring sets for RRT*.

## Design & Implementation
Clamps the radius with `max(Float64(radius_m), 0.0)` and forwards to `hypr_rrt_near_indices(tree, q, radius)`. The returned vector is mutated by the caller (`rpo_rrt_star_add_node!` pushes the nearest index when it is empty), so the delegate must return a fresh `Vector{Int}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RPORRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `q` | SVector{3, Float64} | n/a | yes | Positional argument `q`. |
| in | `radius_m` | Real | n/a | yes | Positional argument `radius_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_near_indices`. Returns `hypr_rrt_near_indices(tree, q, max(Float64(radius_m), 0.0))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_star_add_node_bang|rpo_rrt_star_add_node!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:251-251`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:50-50`
- `callees` → [[gnc.hypr_utils_hypr_rrt_near_indices|hypr_rrt_near_indices]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:50-50`
<!-- vulcan:connections:end -->

## Limitations
A negative radius is silently converted to `0.0`, which returns only coincident nodes rather than signalling a configuration error. The scan is O(N) per call and is invoked every RRT* iteration. Whether the result includes nodes at exactly `radius_m` depends on the delegate's comparison.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 49.
