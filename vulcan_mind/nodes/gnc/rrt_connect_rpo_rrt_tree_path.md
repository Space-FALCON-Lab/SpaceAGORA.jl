---
id: gnc.rrt_connect_rpo_rrt_tree_path
label: rpo_rrt_tree_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_tree_path
  lines:
  - 204
  - 204
inputs:
- id: tree
  type: RPORRTConnectTree
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  description: Return value of `rpo_rrt_tree_path`. Returns `hypr_rrt_tree_path(tree,
    idx)`.
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

# rpo_rrt_tree_path

## Purpose
Reconstructs the ordered sequence of node positions from the tree root to node `idx`, used to extract the raw RRT* solution.

## Design & Implementation
Delegates to `hypr_rrt_tree_path(tree, idx)`, which follows `tree.parents` from `idx` back to the root (parent `0`) and returns the reversed chain as a `3 x k` matrix. `rpo_rrt_star_plan_path` concatenates the goal column onto the result with `hcat`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RPORRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `idx` | Integer | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_tree_path`. Returns `hypr_rrt_tree_path(tree, idx)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:573-573`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_tree_path|hypr_rrt_tree_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:205-205`
<!-- vulcan:connections:end -->

## Limitations
A corrupted `parents` vector containing a cycle would make the delegate loop forever; no cycle guard exists in this layer. `idx` is not bounds-checked before delegation. The returned matrix is freshly allocated on every call, including inside the RRT* improvement loop.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 204.
