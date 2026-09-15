---
id: gnc.planner_core__robot_arm_hypr_post_refine_points
label: _robot_arm_hypr_post_refine_points
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/planner_core.jl
  symbol: _robot_arm_hypr_post_refine_points
  lines:
  - 61
  - 61
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
  description: Return value of `_robot_arm_hypr_post_refine_points`. Returns `current,
    current_components, false, 0` or `current, current_components, improved, rounds_run`.
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

# _robot_arm_hypr_post_refine_points

## Purpose
Coordinate-descent post-processing of the best PSO control-point matrix. For each interior control point and each joint it tries a step in both directions, keeps the change if `_robot_arm_hypr_refinement_better` accepts it, and shrinks the step size each round. Returns the refined `Matrix{Float64}`, its cost components, a Bool `improved` flag and the number of rounds run.

## Design & Implementation
The input is copied into `current` and scored once. If `cfg.refinement_enable` is false, `cfg.refinement_rounds == 0`, or there are no interior points (`size(current, 2) <= 2`), it returns immediately with `(current, components, false, 0)`. Step sizes are `cfg.refinement_step_fraction .* spans` with spans from `joint.upper_rad - joint.lower_rad` clamped to at least `1e-9`. For `round in 1:cfg.refinement_rounds`, columns `2:end-1` and every joint dimension `d` are visited; each candidate is `clamp(current[d,j] + direction*steps[d], lo[d], hi[d])`, skipped if the clamp leaves the value unchanged, and scored with `cost_cutoff` set to the incumbent total. After each round `steps .*= cfg.refinement_shrink`; the loop breaks early when a round makes no change. `current` is mutated in place while iterating.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `obstacles` | AbstractVector{RobotArmSphereObstacle} | n/a | yes | Positional argument `obstacles`. |
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_post_refine_points`. Returns `current, current_components, false, 0` or `current, current_components, improved, rounds_run`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:334-334`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:334-334`

**Downstream**

- `callees` → [[gnc.planner_core__robot_arm_hypr_refinement_better|_robot_arm_hypr_refinement_better]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:99-99`
- `callees` → [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:69-69`
<!-- vulcan:connections:end -->

## Limitations
Complexity is `rounds * (n_points - 2) * n_joints * 2` full cost evaluations, each resampling and running forward kinematics for `cfg.n_samples` configurations. Endpoint columns are never moved, so the start and goal are treated as fixed. Since `_robot_arm_hypr_refinement_better` cannot accept when the incumbent total is `Inf`, refinement silently does nothing on a pruned or infeasible-with-Inf path. There is no `rng` involvement, so the result is deterministic but can stall in a local minimum along the axis-aligned search directions.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/planner_core.jl` line 61.
