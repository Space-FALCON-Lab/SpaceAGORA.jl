---
id: gnc.rrt_warmstart__robot_arm_rrt_segment_is_safe
label: _robot_arm_rrt_segment_is_safe
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_segment_is_safe
  lines:
  - 45
  - 45
inputs:
- id: q_from
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_from`.
- id: q_to
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_to`.
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
  description: Return value of `_robot_arm_rrt_segment_is_safe`. Returns `stats.violation_count
    == 0 && stats.min_clearance + 1.0e-9 >= cfg.safe_distance_`.
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

# _robot_arm_rrt_segment_is_safe

## Purpose
`_robot_arm_rrt_segment_is_safe` decides whether a candidate joint-space edge keeps the whole arm clear of the sphere obstacles by the configured safe distance. It gates every extension in `_robot_arm_rrt_extend!`, every shortcut in `_robot_arm_rrt_shortcut_path`, and the initial direct-path test.

## Design & Implementation
Signature `(q_from, q_to, model::ClothArmModel, base_pose::ClothArmBasePose, obstacles::AbstractVector{RobotArmSphereObstacle}, cfg::RobotArmHYPRConfig)`. With no obstacles it returns `true` immediately. Otherwise it builds samples with `_robot_arm_rrt_segment_samples(q_from, q_to, cfg.rrt_warmstart_collision_sample_ds_rad)`, evaluates `robot_arm_clearance_stats_from_samples(model, base_pose, samples, obstacles, cfg.safe_distance_m)` (forward kinematics per sample and minimum link-to-sphere distance), and returns `stats.violation_count == 0 && stats.min_clearance + 1e-9 >= cfg.safe_distance_m`. The `1e-9` slack prevents rejecting edges that touch the safe distance exactly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_from` | Any | n/a | yes | Positional argument `q_from`. |
| in | `q_to` | Any | n/a | yes | Positional argument `q_to`. |
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `obstacles` | AbstractVector{RobotArmSphereObstacle} | n/a | yes | Positional argument `obstacles`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_segment_is_safe`. Returns `stats.violation_count == 0 && stats.min_clearance + 1.0e-9 >= cfg.safe_distance_`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_warmstart__robot_arm_rrt_extend_bang|_robot_arm_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:71-71`
- [[gnc.rrt_warmstart__robot_arm_rrt_shortcut_path|_robot_arm_rrt_shortcut_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:125-125`
- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:229-229`

**Downstream**

- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_segment_samples|_robot_arm_rrt_segment_samples]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:54-54`
- `callees` → [[gncz.clearance_robot_arm_clearance_stats_from_samples|robot_arm_clearance_stats_from_samples]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:55-55`
<!-- vulcan:connections:end -->

## Limitations
Safety is verified only at discrete samples; the arm can sweep through an obstacle between samples if `rrt_warmstart_collision_sample_ds_rad` is coarse relative to obstacle size. Self-collision and joint limits are not checked here. The cost is a full forward-kinematics evaluation per sample, which makes this the dominant expense of the warm start. Only spherical obstacles are supported.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 45.
