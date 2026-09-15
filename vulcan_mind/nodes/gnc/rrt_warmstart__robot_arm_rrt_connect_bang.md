---
id: gnc.rrt_warmstart__robot_arm_rrt_connect_bang
label: _robot_arm_rrt_connect!
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_connect!
  lines:
  - 81
  - 81
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
  description: Return value of `_robot_arm_rrt_connect!`; mutates `tree` in place.
    Returns `status, idx`.
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

# _robot_arm_rrt_connect!

## Purpose
`_robot_arm_rrt_connect!` greedily extends a tree toward a fixed target configuration step after step until the target is reached, the extension is blocked, or a step budget is exhausted. It is the Connect half of RRT-Connect, applied to the opposite tree after each exploration extension.

## Design & Implementation
Signature `(tree::RobotArmRRTConnectTree, q_target::AbstractVector{<:Real}, model, base_pose, obstacles, cfg::RobotArmHYPRConfig)`. It initialises `status = :advanced`, `idx = 1`, `steps = 0` and loops `while status == :advanced && steps < cfg.rrt_warmstart_connect_max_steps`, calling `_robot_arm_rrt_extend!` and incrementing `steps`. It returns the final `(status, idx)`, where `status` is `:reached` when the target was attained, `:trapped` when blocked, or `:advanced` when the step budget ran out. The tree is mutated through the extension calls.

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
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_connect!`; mutates `tree` in place. Returns `status, idx`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:265-265`

**Downstream**

- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_extend_bang|_robot_arm_rrt_extend!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:93-93`
<!-- vulcan:connections:end -->

## Limitations
Hitting `rrt_warmstart_connect_max_steps` returns `:advanced`, which the caller treats as not connected; the partial branch stays in the tree and inflates future nearest-neighbour scans. Every step repeats a linear nearest search even though the nearest node is almost always the one just added. If `rrt_warmstart_connect_max_steps` is zero the loop never runs and `idx = 1` (the root) is returned with `:advanced`.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 81.
