---
id: gnc.rrt_connect_rpo_rrt_connect_join_paths
label: rpo_rrt_connect_join_paths
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_connect_join_paths
  lines:
  - 214
  - 214
inputs:
- id: start_tree
  type: RPORRTConnectTree
  units: n/a
  required: true
  description: Positional argument `start_tree`.
- id: start_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `start_idx`.
- id: goal_tree
  type: RPORRTConnectTree
  units: n/a
  required: true
  description: Positional argument `goal_tree`.
- id: goal_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `goal_idx`.
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
  description: Return value of `rpo_rrt_connect_join_paths`. Returns `hypr_rrt_join_paths(start_tree,
    start_idx, goal_tree, goal_idx)`.
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

# rpo_rrt_connect_join_paths

## Purpose
Stitches the start-tree path and the goal-tree path together at the connection point into one start-to-goal polyline for RRT-Connect.

## Design & Implementation
Delegates to `hypr_rrt_join_paths(start_tree, start_idx, goal_tree, goal_idx)`, which traces `start_idx` to the start root, traces `goal_idx` to the goal root, reverses the latter, and concatenates. `rpo_rrt_connect_plan_path` selects which tree index is the start side based on the parity of the iteration that achieved `:reached`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_tree` | RPORRTConnectTree | n/a | yes | Positional argument `start_tree`. |
| in | `start_idx` | Integer | n/a | yes | Positional argument `start_idx`. |
| in | `goal_tree` | RPORRTConnectTree | n/a | yes | Positional argument `goal_tree`. |
| in | `goal_idx` | Integer | n/a | yes | Positional argument `goal_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_connect_join_paths`. Returns `hypr_rrt_join_paths(start_tree, start_idx, goal_tree, goal_idx)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:390-390`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_join_paths|hypr_rrt_join_paths]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:215-215`
<!-- vulcan:connections:end -->

## Limitations
The junction usually contains two nearly coincident points (the last node of each tree), producing a zero-length segment that later shortcutting must remove. Index validity is assumed from the caller. The output is a dense matrix built from two tree traversals and is allocated once per successful plan.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 214.
