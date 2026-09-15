---
id: gncz.planner_core_plan_robot_arm_motion_hypr
label: plan_robot_arm_motion_hypr
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/planner_core.jl
  symbol: plan_robot_arm_motion_hypr
  lines:
  - 174
  - 364
inputs:
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
- id: q_start
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_start`.
- id: target
  type: Any
  units: n/a
  required: true
  description: Positional argument `target`.
- id: planner_config
  type: RobotArmPlannerConfig
  units: n/a
  required: false
  description: Keyword argument `planner_config` (default `RobotArmPlannerConfig()`).
- id: hypr_config
  type: RobotArmHYPRConfig
  units: n/a
  required: false
  description: Keyword argument `hypr_config` (default `RobotArmHYPRConfig()`).
- id: obstacles
  type: AbstractVector{RobotArmSphereObstacle}
  units: n/a
  required: false
  description: Keyword argument `obstacles` (default `RobotArmSphereObstacle[]`).
- id: rng
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rng` (default `Random.default_rng()`).
- id: module_api
  type: Module
  units: n/a
  required: true
  description: RobotArmPlanning namespace providing inverse kinematics, the configuration
    types, and the swarm, clearance, and warmstart routines.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: RobotArmHYPRResult
  units: n/a
  description: Return value of `plan_robot_arm_motion_hypr`. Returns `RobotArmHYPRResult(plan,
    points, sampled, comps.total, comps, [comps.total], cfg` or `robot_arm_hypr_path_cost_components(points,
    model, base_pose, obstacles, cfg; co` or `RobotArmHYPRResult(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# plan_robot_arm_motion_hypr

## Purpose
`plan_robot_arm_motion_hypr` is the entry point of the hybrid robot-arm planner. It takes a start configuration and a Cartesian target and returns a complete joint-space plan together with the diagnostics needed to judge whether the plan is trustworthy.

## Theory & Math
The search optimises the interior control points of a curve in joint space. A particle is a flattened vector of $n \times w$ joint values for $w$ interior waypoints, bounded by the joint limits repeated per waypoint. The cost combines path length, smoothness, and the squared clearance deficit, and the swarm updates each particle with the usual inertia and cognitive and social attraction terms, with the coefficients annealed by the configured schedule.

## Model & Assumptions
The goal configuration comes from damped least-squares inverse kinematics seeded with the start configuration, so the plan is only as reachable as that solve. A degenerate configuration with zero interior waypoints short circuits to a straight polyline between start and goal, still retimed and scored so the returned record has the same shape as a full run.

## Design & Implementation
An optional bidirectional rapidly exploring random tree supplies a feasible seed path, resampled to the waypoint count; without it the swarm is seeded by linear interpolation. The first particle takes the seed exactly and the rest are scattered by a spread fraction of the joint span. After the swarm converges the best point set is locally refined, retimed under velocity, acceleration, and base wrench limits, and converted into a plan whose components record refinement, retiming, and warmstart outcomes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ClothArmModel | n/a | yes | Positional argument `model`. |
| in | `base_pose` | ClothArmBasePose | n/a | yes | Positional argument `base_pose`. |
| in | `q_start` | Any | n/a | yes | Positional argument `q_start`. |
| in | `target` | Any | n/a | yes | Positional argument `target`. |
| in | `planner_config` | RobotArmPlannerConfig | n/a | no | Keyword argument `planner_config` (default `RobotArmPlannerConfig()`). |
| in | `hypr_config` | RobotArmHYPRConfig | n/a | no | Keyword argument `hypr_config` (default `RobotArmHYPRConfig()`). |
| in | `obstacles` | AbstractVector{RobotArmSphereObstacle} | n/a | no | Keyword argument `obstacles` (default `RobotArmSphereObstacle[]`). |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `Random.default_rng()`). |
| in | `module_api` | Module | n/a | yes | RobotArmPlanning namespace providing inverse kinematics, the configuration types, and the swarm, clearance, and warmstart routines. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmHYPRResult | n/a | — | Return value of `plan_robot_arm_motion_hypr`. Returns `RobotArmHYPRResult(plan, points, sampled, comps.total, comps, [comps.total], cfg` or `robot_arm_hypr_path_cost_components(points, model, base_pose, obstacles, cfg; co` or `RobotArmHYPRResult(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.robot_arm_planning_plan_robot_arm_motion|plan_robot_arm_motion]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:84-84`

**Downstream**

- `callees` → [[gnc.config__validate_robot_arm_hypr_config|_validate_robot_arm_hypr_config]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:184-184`
- `callees` → [[gnc.config_robotarmhyprresult|RobotArmHYPRResult]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:218-218`
- `callees` → [[gnc.planner_core__robot_arm_hypr_post_refine_points|_robot_arm_hypr_post_refine_points]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:334-334`
- `callees` → [[gnc.planner_core__robot_arm_plan_from_q_reference|_robot_arm_plan_from_q_reference]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:204-204`
- `callees` → [[gnc.planner_core_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:268-268`
- `callees` → [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:206-206`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:289-289`
- `callees` → [[gnc.pso_path_planning_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:268-268`
- `callees` → [[gnc.robot_arm_planning__reference_times|_reference_times]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:201-201`
- `callees` → [[gnc.robot_arm_planning_robotarmplannerconfig|RobotArmPlannerConfig]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:179-179`
- `callees` → [[gnc.rrt_warmstart__robot_arm_empty_rrt_warmstart_diagnostics|_robot_arm_empty_rrt_warmstart_diagnostics]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:197-197`
- `callees` → [[gnc.rrt_warmstart__robot_arm_resample_polyline_points|_robot_arm_resample_polyline_points]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:242-242`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_warmstart_fields|_robot_arm_rrt_warmstart_fields]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:217-217`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_control_points|_robot_arm_control_points]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:269-269`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_flatten_internal_points|_robot_arm_flatten_internal_points]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:245-245`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_base_wrench_ratios|_robot_arm_hypr_base_wrench_ratios]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:205-205`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_cull_swarm_bang|_robot_arm_hypr_cull_swarm!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:306-306`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_early_stopping_feasible|_robot_arm_hypr_early_stopping_feasible]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:291-291`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_iteration_weights|_robot_arm_hypr_iteration_weights]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:319-319`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_material_improvement|_robot_arm_hypr_material_improvement]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:292-292`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_seed_control_points|_robot_arm_seed_control_points]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:241-241`
- `callees` → [[gnc.swarm_and_retiming_robot_arm_sample_hypr_path|robot_arm_sample_hypr_path]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:200-200`
- `callees` → [[gncz.config_robotarmhyprconfig|RobotArmHYPRConfig]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:180-180`
- `callees` → [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:239-239`
- `callees` → [[gncz.swarm_and_retiming__robot_arm_hypr_retime_reference|_robot_arm_hypr_retime_reference]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:203-203`
- `callees` → [[vehicle.robotics_cloth_ik|cloth_ik]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:188-188`
<!-- vulcan:connections:end -->

## Limitations
The optimiser is stochastic, so results depend on the supplied random generator, and neither the swarm nor the refinement guarantees a collision-free path. Cost evaluation is the bottleneck because it runs forward kinematics at every sample of every particle at every iteration.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/planner_core.jl:1-364`.
