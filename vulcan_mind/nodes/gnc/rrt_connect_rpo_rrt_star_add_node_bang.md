---
id: gnc.rrt_connect_rpo_rrt_star_add_node_bang
label: rpo_rrt_star_add_node!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_star_add_node!
  lines:
  - 243
  - 243
inputs:
- id: tree
  type: RPORRTConnectTree
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: q_new
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `q_new`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: settings
  type: RPORRTStarSettings
  units: n/a
  required: true
  description: Positional argument `settings`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `safe_distance_m`.
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
  type: Int
  units: n/a
  description: Return value of `rpo_rrt_star_add_node!`; mutates `tree` in place.
    Returns `0` or `new_idx`.
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

# rpo_rrt_star_add_node!

## Purpose
Inserts `q_new` into an RRT* tree by choosing the cheapest collision-free parent among nearby nodes, then rewires neighbours through the new node when that lowers their cost.

## Design & Implementation
Gets `nearest_idx` and `near_idxs` (radius `settings.neighbor_radius_m`, seeded with the nearest index when empty). Initial `best_cost` is via the nearest node; each near candidate is adopted when `costs[idx] + edge_cost + 1e-9 < best_cost` and its edge passes `rpo_rrt_segment_is_safe`. The chosen parent edge is collision-checked again; on failure the function returns `0` without mutating the tree. Otherwise it `push!`es the node, parent, and cost, then for every other near node tests `best_cost + norm(nodes[idx] - q_new) + 1e-9 < costs[idx]` plus edge safety, and on success sets `parents[idx] = new_idx`, `costs[idx] = candidate_cost`, and calls `rpo_rrt_refresh_subtree_costs!(tree, idx)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RPORRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `q_new` | SVector{3, Float64} | n/a | yes | Positional argument `q_new`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `settings` | RPORRTStarSettings | n/a | yes | Positional argument `settings`. |
| in | `safe_distance_m` | Real | n/a | yes | Keyword argument `safe_distance_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `rpo_rrt_star_add_node!`; mutates `tree` in place. Returns `0` or `new_idx`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:552-552`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:252-252`
- `callees` → [[gnc.rrt_connect_rpo_rrt_near_indices|rpo_rrt_near_indices]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:251-251`
- `callees` → [[gnc.rrt_connect_rpo_rrt_nearest_index|rpo_rrt_nearest_index]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:250-250`
- `callees` → [[gnc.rrt_connect_rpo_rrt_refresh_subtree_costs_bang|rpo_rrt_refresh_subtree_costs!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:301-301`
- `callees` → [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:260-260`
<!-- vulcan:connections:end -->

## Limitations
The nearest-node edge is checked twice when it wins (once in the candidate loop only if it appears in `near_idxs`, then again unconditionally), wasting a segment check. Return value `0` is a sentinel rather than a `Nothing`, so callers must compare explicitly. Rewiring uses the Euclidean edge length only, ignoring the clearance and smoothness terms that `rpo_normalized_path_cost_components` later scores.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 243.
