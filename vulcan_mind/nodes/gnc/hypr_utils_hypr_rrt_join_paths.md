---
id: gnc.hypr_utils_hypr_rrt_join_paths
label: hypr_rrt_join_paths
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_rrt_join_paths
  lines:
  - 189
  - 189
inputs:
- id: start_tree
  type: Any
  units: n/a
  required: true
  description: Positional argument `start_tree`.
- id: start_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `start_idx`.
- id: goal_tree
  type: Any
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
  description: Return value of `hypr_rrt_join_paths`. Returns `hcat(start_path, reverse(goal_path[:,
    1:(end - 1)]; dims=2))`.
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

# hypr_rrt_join_paths

## Purpose
Splices the start-tree path and the goal-tree path together at the point where RRT-Connect joined them, producing one start-to-goal state sequence.

## Design & Implementation
Reconstructs both root-to-node paths, then concatenates the start path with the goal path reversed along its columns after dropping the goal path's last column, which duplicates the connection state already present at the end of the start path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_tree` | Any | n/a | yes | Positional argument `start_tree`. |
| in | `start_idx` | Integer | n/a | yes | Positional argument `start_idx`. |
| in | `goal_tree` | Any | n/a | yes | Positional argument `goal_tree`. |
| in | `goal_idx` | Integer | n/a | yes | Positional argument `goal_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_rrt_join_paths`. Returns `hcat(start_path, reverse(goal_path[:, 1:(end - 1)]; dims=2))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_connect_join_paths|rpo_rrt_connect_join_paths]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:215-215`
- [[gnc.rrt_warmstart__robot_arm_rrt_join_paths|_robot_arm_rrt_join_paths]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:106-106`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_tree_path|hypr_rrt_tree_path]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:190-190`
<!-- vulcan:connections:end -->

## Limitations
Assumes the two nodes are the same state, or within the connect tolerance; if the trees met with a residual gap the joined path has a small discontinuity at the seam that later shortcutting is expected to smooth.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 189.
