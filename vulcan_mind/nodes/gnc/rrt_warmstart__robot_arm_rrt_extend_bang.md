---
id: gnc.rrt_warmstart__robot_arm_rrt_extend_bang
label: _robot_arm_rrt_extend!
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_extend!
  lines:
  - 60
  - 60
inputs:
- id: tree
  type: RobotArmRRTConnectTree
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: q_target
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `q_target`.
- id: model
  type: ClothArmModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: base_pose
  type: ClothArmBasePose
  units: n/a
  required: true
  description: Positional argument `base_pose`.
- id: obstacles
  type: AbstractVector{RobotArmSphereObstacle}
  units: n/a
  required: true
  description: Positional argument `obstacles`.
- id: cfg
  type: RobotArmHYPRConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_robot_arm_rrt_extend!`; mutates `tree` in place.
    Returns `:trapped, nearest_idx` or `status, length(tree.nodes)`.
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

# _robot_arm_rrt_extend!

## Purpose
`_robot_arm_rrt_extend!` performs one RRT extension: find the nearest tree node to a target, steer one step toward it, verify the new edge is collision-free, and append the new node. It is the primitive both the exploration step and `_robot_arm_rrt_connect!` are built on.

## Design & Implementation
Signature `(tree::RobotArmRRTConnectTree, q_target::AbstractVector{<:Real}, model, base_pose, obstacles, cfg::RobotArmHYPRConfig)`. It calls `_robot_arm_rrt_nearest_index`, then `_robot_arm_rrt_steer(tree.nodes[nearest_idx], q_target, cfg.rrt_warmstart_step_size_rad)`. If steering reports `:trapped`, or `_robot_arm_rrt_segment_is_safe` rejects the edge, it returns `(:trapped, nearest_idx)` without modifying the tree. Otherwise it pushes `q_new` onto `tree.nodes`, `nearest_idx` onto `tree.parents`, and `tree.costs[nearest_idx] + norm(q_new - nodes[nearest_idx])` onto `tree.costs`, returning `(status, length(tree.nodes))` where `status` is `:advanced` or `:reached`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RobotArmRRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `q_target` | AbstractVector{<:Real} | n/a | yes | Positional argument `q_target`. |
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `obstacles` | AbstractVector{RobotArmSphereObstacle} | n/a | yes | Positional argument `obstacles`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_extend!`; mutates `tree` in place. Returns `:trapped, nearest_idx` or `status, length(tree.nodes)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_warmstart__robot_arm_rrt_connect_bang|_robot_arm_rrt_connect!]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:93-93`
- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:261-261`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:74-74`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_nearest_index|_robot_arm_rrt_nearest_index]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:68-68`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_segment_is_safe|_robot_arm_rrt_segment_is_safe]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:71-71`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_steer|_robot_arm_rrt_steer]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:69-69`
<!-- vulcan:connections:end -->

## Limitations
A `:reached` status means the steered point equals `q_target`, but the returned index is the new node's, so callers must treat the node as the target. The tree grows without bound; there is no node cap other than the caller's iteration limit. Each call performs a linear nearest-neighbour scan and a full sampled collision check. The three tree vectors are mutated non-atomically, so an exception between the pushes would leave them inconsistent (none of the called functions is expected to throw after the first push).

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 60.
