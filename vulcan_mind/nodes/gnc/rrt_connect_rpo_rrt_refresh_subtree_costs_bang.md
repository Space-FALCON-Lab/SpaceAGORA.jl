---
id: gnc.rrt_connect_rpo_rrt_refresh_subtree_costs_bang
label: rpo_rrt_refresh_subtree_costs!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_refresh_subtree_costs!
  lines:
  - 209
  - 209
inputs:
- id: tree
  type: RPORRTConnectTree
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: parent_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `parent_idx`.
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
  description: Return value of `rpo_rrt_refresh_subtree_costs!`; mutates `tree` in
    place. Returns `hypr_rrt_refresh_subtree_costs!(tree, parent_idx)`.
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

# rpo_rrt_refresh_subtree_costs!

## Purpose
After RRT* rewires a node to a cheaper parent, propagates the updated accumulated cost to every descendant of `parent_idx`.

## Design & Implementation
Forwards to `hypr_rrt_refresh_subtree_costs!(tree, parent_idx)`, which mutates `tree.costs` in place for all nodes whose ancestry passes through `parent_idx`. Invoked from `rpo_rrt_star_add_node!` once per rewired neighbour so that later cost comparisons see consistent values.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RPORRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `parent_idx` | Integer | n/a | yes | Positional argument `parent_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_refresh_subtree_costs!`; mutates `tree` in place. Returns `hypr_rrt_refresh_subtree_costs!(tree, parent_idx)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_star_add_node_bang|rpo_rrt_star_add_node!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:301-301`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_refresh_subtree_costs_bang|hypr_rrt_refresh_subtree_costs!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:210-210`
<!-- vulcan:connections:end -->

## Limitations
The subtree walk is linear in tree size for each rewire, so a rewiring-heavy iteration can cost O(N * near-set size). Only `costs` is refreshed; `parents` must already be correct. The function assumes the parent's own cost was updated by the caller before the call.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 209.
