---
id: gnc.rrt_warmstart__robot_arm_rrt_shortcut_path
label: _robot_arm_rrt_shortcut_path
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_shortcut_path
  lines:
  - 110
  - 110
inputs:
- id: path
  type: Any
  units: n/a
  required: true
  description: Positional argument `path`.
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
- id: rng
  type: Any
  units: n/a
  required: true
  description: Positional argument `rng`.
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
  description: Return value of `_robot_arm_rrt_shortcut_path`. Returns `pts`.
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

# _robot_arm_rrt_shortcut_path

## Purpose
`_robot_arm_rrt_shortcut_path` post-processes a jagged RRT path by repeatedly trying to replace a random sub-sequence of waypoints with a single straight joint-space segment, keeping the replacement only when it remains collision-free and does not lengthen the path. It shortens the warm start before it is scored and handed to HYPR.

## Design & Implementation
Signature `(path, model, base_pose, obstacles, cfg::RobotArmHYPRConfig, rng)`. The path is copied to `pts = Matrix{Float64}(path)`; paths with two or fewer columns are returned as-is. For `cfg.rrt_warmstart_shortcut_iters` iterations it draws `i in 1:(n_pts-2)` and `j in (i+2):n_pts` (so at least one interior waypoint is removed), checks `_robot_arm_rrt_segment_is_safe(pts[:, i], pts[:, j], ...)`, forms `candidate = pts[:, vcat(1:i, j:n_pts)]`, and adopts it when `_robot_arm_path_length(candidate) <= _robot_arm_path_length(pts) + 1e-9`. The loop breaks early once only two points remain. Returns the reduced matrix.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | Any | n/a | yes | Positional argument `path`. |
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `obstacles` | AbstractVector{RobotArmSphereObstacle} | n/a | yes | Positional argument `obstacles`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `rng` | Any | n/a | yes | Positional argument `rng`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_shortcut_path`. Returns `pts`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:277-277`

**Downstream**

- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_segment_is_safe|_robot_arm_rrt_segment_is_safe]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:125-125`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_path_length|_robot_arm_path_length]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:127-127`
<!-- vulcan:connections:end -->

## Limitations
Random pair selection means the outcome depends on `rng` and the fixed iteration count; no convergence criterion exists, and with few waypoints many iterations are wasted on repeats. The straight-line replacement is checked at sample resolution only, inheriting the discretisation gaps of `_robot_arm_rrt_segment_is_safe`. Because a straight segment between two points is never longer than the polyline it replaces, the length test is effectively always true; it exists only to guard numerical noise. Each iteration allocates a new candidate matrix.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 110.
