---
id: gnc.rrt_warmstart__robot_arm_rrt_path_score
label: _robot_arm_rrt_path_score
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_path_score
  lines:
  - 169
  - 169
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
  description: Return value of `_robot_arm_rrt_path_score`. Returns `(_robot_arm_path_length(samples)
    / len_ref)^2 +`.
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

# _robot_arm_rrt_path_score

## Purpose
`_robot_arm_rrt_path_score` evaluates a candidate warm-start path with the same objective components HYPR optimises (normalised path length squared plus weighted obstacle violation and clearance penalty), so the RRT result and the optimiser's later iterates are comparable and the diagnostics report a meaningful `cost`.

## Theory & Math
$$J = \Big(\frac{L_{samples}}{\max(\lVert q_f - q_0\rVert, 10^{-6})}\Big)^2 + w_{obs}\Big(N_{viol} + \frac{P_{clear}}{\max(d_{safe}^2, 10^{-8})}\Big)$$

where $L_{samples}$ is the resampled joint-space path length (rad), $q_0, q_f$ the endpoints, $w_{obs}$ = `cfg.w_obs`, $N_{viol}$ the number of samples violating clearance, $P_{clear}$ the accumulated clearance penalty (m^2) and $d_{safe}$ = `cfg.safe_distance_m`.

## Design & Implementation
Signature `(path, model, base_pose, obstacles, cfg::RobotArmHYPRConfig)`. It resamples the path with `robot_arm_sample_hypr_path(path, cfg.n_samples; curve_type=:polyline)`, sets `len_ref = max(norm(path[:, end] - path[:, 1]), 1e-6)` (straight-line joint distance), computes clearance statistics via `robot_arm_clearance_stats_from_samples(..., cfg.safe_distance_m)`, and with `penalty_ref = max(cfg.safe_distance_m^2, 1e-8)` returns `(path_length(samples) / len_ref)^2 + cfg.w_obs * (violation_count + clearance_penalty / penalty_ref)`. A perfectly straight collision-free path therefore scores exactly 1.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | Any | n/a | yes | Positional argument `path`. |
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `obstacles` | AbstractVector{RobotArmSphereObstacle} | n/a | yes | Positional argument `obstacles`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_path_score`. Returns `(_robot_arm_path_length(samples) / len_ref)^2 +`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:230-230`

**Downstream**

- `callees` → [[gnc.swarm_and_retiming__robot_arm_path_length|_robot_arm_path_length]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:180-180`
- `callees` → [[gnc.swarm_and_retiming_robot_arm_sample_hypr_path|robot_arm_sample_hypr_path]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:176-176`
- `callees` → [[gncz.clearance_robot_arm_clearance_stats_from_samples|robot_arm_clearance_stats_from_samples]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:178-178`
<!-- vulcan:connections:end -->

## Limitations
The score omits any HYPR terms other than length and obstacle clearance (for example smoothness or joint-effort weights), so it is only a partial proxy for the optimiser's objective. `len_ref` collapses to `1e-6` when start and goal coincide, making the length term blow up for any non-trivial path. The path is resampled to `cfg.n_samples` before measuring, so very short polylines with sharp corners are measured on the resampled approximation.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 169.
