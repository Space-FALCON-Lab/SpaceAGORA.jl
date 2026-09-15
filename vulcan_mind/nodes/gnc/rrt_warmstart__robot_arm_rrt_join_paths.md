---
id: gnc.rrt_warmstart__robot_arm_rrt_join_paths
label: _robot_arm_rrt_join_paths
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_join_paths
  lines:
  - 100
  - 100
inputs:
- id: start_tree
  type: RobotArmRRTConnectTree
  units: n/a
  required: true
  description: Positional argument `start_tree`.
- id: start_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `start_idx`.
- id: goal_tree
  type: RobotArmRRTConnectTree
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
  description: Return value of `_robot_arm_rrt_join_paths`. Returns `hypr_rrt_join_paths(start_tree,
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

# _robot_arm_rrt_join_paths

## Purpose
`_robot_arm_rrt_join_paths` assembles the final joint-space path once the start and goal trees have met, by concatenating the root-to-meeting-node path of the start tree with the reversed root-to-meeting-node path of the goal tree.

## Design & Implementation
Signature `(start_tree::RobotArmRRTConnectTree, start_idx::Integer, goal_tree::RobotArmRRTConnectTree, goal_idx::Integer)`. It forwards to `hypr_rrt_join_paths`, which calls `hypr_rrt_tree_path` on each tree (walking `parents` from the index back to the root and returning a matrix with one configuration per column, root first) and returns `hcat(start_path, reverse(goal_path[:, 1:end-1]; dims=2))`. Dropping the last column of `goal_path` before reversing removes the duplicated meeting configuration. The output is a `Matrix{Float64}` of size `n_joints × n_waypoints` running from `q_start` to `q_goal`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_tree` | RobotArmRRTConnectTree | n/a | yes | Positional argument `start_tree`. |
| in | `start_idx` | Integer | n/a | yes | Positional argument `start_idx`. |
| in | `goal_tree` | RobotArmRRTConnectTree | n/a | yes | Positional argument `goal_tree`. |
| in | `goal_idx` | Integer | n/a | yes | Positional argument `goal_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_join_paths`. Returns `hypr_rrt_join_paths(start_tree, start_idx, goal_tree, goal_idx)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:269-269`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_join_paths|hypr_rrt_join_paths]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:106-106`
<!-- vulcan:connections:end -->

## Limitations
Correctness assumes `start_tree.nodes[start_idx]` and `goal_tree.nodes[goal_idx]` are the same configuration; the function does not verify that, so a caller passing mismatched indices gets a path with a hidden jump. Path reconstruction walks the parent chain, which is O(depth) with a vector growing per step. The result is unsimplified and typically jagged, hence the subsequent shortcut pass.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 100.
