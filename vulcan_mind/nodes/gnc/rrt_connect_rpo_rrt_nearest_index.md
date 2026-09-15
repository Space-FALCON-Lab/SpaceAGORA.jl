---
id: gnc.rrt_connect_rpo_rrt_nearest_index
label: rpo_rrt_nearest_index
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_nearest_index
  lines:
  - 44
  - 44
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
  description: Return value of `rpo_rrt_nearest_index`. Returns `hypr_rrt_nearest_index(tree,
    q)`.
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

# rpo_rrt_nearest_index

## Purpose
Returns the index of the tree node closest (Euclidean) to a query point `q`, used as the branch point for every extend step.

## Design & Implementation
A thin typed wrapper: `rpo_rrt_nearest_index(tree::RPORRTConnectTree, q::SVector{3,Float64})` forwards to the shared `hypr_rrt_nearest_index(tree, q)`, which performs a brute-force scan over `tree.nodes`. Called once per extension by `rpo_rrt_extend!`, `rpo_rrt_star_add_node!`, and directly by `rpo_rrt_star_plan_path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RPORRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `q` | SVector{3, Float64} | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_nearest_index`. Returns `hypr_rrt_nearest_index(tree, q)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_extend_bang|rpo_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:161-161`
- [[gnc.rrt_connect_rpo_rrt_star_add_node_bang|rpo_rrt_star_add_node!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:250-250`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:547-547`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_nearest_index|hypr_rrt_nearest_index]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:45-45`
<!-- vulcan:connections:end -->

## Limitations
Linear-time search with no spatial index, making each planner iteration O(N) in tree size. Ties are resolved by whatever order the delegate uses (first minimum). An empty tree is impossible via the public constructor but would make the delegate return an invalid index.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 44.
