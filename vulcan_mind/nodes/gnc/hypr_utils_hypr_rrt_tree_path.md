---
id: gnc.hypr_utils_hypr_rrt_tree_path
label: hypr_rrt_tree_path
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_rrt_tree_path
  lines:
  - 176
  - 176
inputs:
- id: tree
  type: Any
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
  description: Return value of `hypr_rrt_tree_path`. Returns `reduce(hcat, nodes)`.
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

# hypr_rrt_tree_path

## Purpose
Reconstructs the sequence of states from the tree root to a given node by following parent links.

## Design & Implementation
Starts an empty vector typed to the node element type, walks `current = parents[current]` until a non-positive index marks the root, pushing each node, then reverses in place and concatenates the states into a matrix with `reduce(hcat, ...)`. The parent accessor handles either tree layout.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | Any | n/a | yes | Positional argument `tree`. |
| in | `idx` | Integer | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_rrt_tree_path`. Returns `reduce(hcat, nodes)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.hypr_utils_hypr_rrt_join_paths|hypr_rrt_join_paths]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:190-190`
- [[gnc.rrt_connect_rpo_rrt_tree_path|rpo_rrt_tree_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:205-205`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils__hypr_rrt_parents|_hypr_rrt_parents]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:179-179`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:181-181`
<!-- vulcan:connections:end -->

## Limitations
A cycle in the parent vector, which a buggy rewiring could create, makes this loop forever; there is no visited set or depth cap.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 176.
