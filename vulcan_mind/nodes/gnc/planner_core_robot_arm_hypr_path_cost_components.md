---
id: gnc.planner_core_robot_arm_hypr_path_cost_components
label: robot_arm_hypr_path_cost_components
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/planner_core.jl
  symbol: robot_arm_hypr_path_cost_components
  lines:
  - 2
  - 2
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
- id: cost_cutoff
  type: Real
  units: n/a
  required: false
  description: Keyword argument `cost_cutoff` (default `Inf`).
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
  description: Return value of `robot_arm_hypr_path_cost_components`. Returns `(`.
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

# robot_arm_hypr_path_cost_components

## Purpose
Scores a candidate robot-arm joint path, given as a matrix `points` (joint index by control point, radians), against the HYPR cost model. It returns a named tuple with the weighted `total` plus every component (`J_len`, `J_len_norm`, `J_smooth`, `J_obs`, clearance statistics) so the PSO search, refinement loop and result diagnostics can share one evaluation.

## Theory & Math
$$J = w_{len}\left(\frac{L}{L_{ref}}\right)^2 + w_{smooth} J_{smooth} + w_{obs}\left(N_{viol} + \frac{P_{clr}}{\max(d_{safe}^2, 10^{-8})}\right)$$ where $L$ is the sampled path length in joint space, $L_{ref} = \max(\|q_{end} - q_{start}\|, 10^{-6})$, $N_{viol}$ is the number of samples violating clearance, $P_{clr}$ the accumulated clearance penalty, and $d_{safe}$ = `cfg.safe_distance_m`.

## Design & Implementation
The control points are resampled into `cfg.n_samples` joint configurations with `robot_arm_sample_hypr_path` using `cfg.curve_type`. Length is normalised by `len_ref = max(norm(points[:, end] - points[:, 1]), 1e-6)` so `J_len_norm` is scale-free; smoothness comes from `_robot_arm_path_smoothness`. Obstacle cost uses `robot_arm_clearance_stats_from_samples` against `obstacles` (spheres) with `cfg.safe_distance_m`, forming `J_obs = violation_count + clearance_penalty / max(safe_distance_m^2, 1e-8)`. The weighted sum is `w_len * J_len_norm^2 + w_smooth * J_smooth + w_obs * J_obs`. If the sum exceeds `cost_cutoff` the tuple is returned with `total = Inf` but all components intact, which lets callers prune particles that cannot beat their personal best.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `obstacles` | AbstractVector{RobotArmSphereObstacle} | n/a | yes | Positional argument `obstacles`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `cost_cutoff` | Real | n/a | no | Keyword argument `cost_cutoff` (default `Inf`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `robot_arm_hypr_path_cost_components`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core__robot_arm_hypr_post_refine_points|_robot_arm_hypr_post_refine_points]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:69-69`
- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:270-270`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:206-206`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:19-19`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_path_length|_robot_arm_path_length]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:12-12`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_path_smoothness|_robot_arm_path_smoothness]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:14-14`
- `callees` → [[gnc.swarm_and_retiming_robot_arm_sample_hypr_path|robot_arm_sample_hypr_path]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:10-10`
- `callees` → [[gncz.clearance_robot_arm_clearance_stats_from_samples|robot_arm_clearance_stats_from_samples]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
The cutoff is applied only after all components have been computed, so `cost_cutoff` saves no forward-kinematics work; it only changes the reported `total`. Path length is measured in raw joint-angle units, mixing joints of different reach, and `len_ref` collapses to `1e-6` when start and goal coincide, inflating `J_len_norm`. Obstacles must be `RobotArmSphereObstacle`; no other geometry is accepted. No validation is done that `points` has one column per joint or at least two columns.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/planner_core.jl` line 2.
