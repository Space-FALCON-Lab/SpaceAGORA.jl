---
id: gnc.planner_core_evaluate_position
label: evaluate_position
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/planner_core.jl
  symbol: evaluate_position
  lines:
  - 268
  - 268
inputs:
- id: pos
  type: Any
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: cost_cutoff
  type: Any
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
  description: Return value of `evaluate_position`. Returns `robot_arm_hypr_path_cost_components(points,
    model, base_pose, obstacles, cfg; co`.
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

# evaluate_position

## Purpose
Closure inside `plan_robot_arm_motion_hypr` that turns one flattened PSO particle `pos` (a vector of length `n_joints * cfg.n_waypoints`, radians) into a control-point matrix and returns its HYPR cost components. It is the single objective function used in every PSO iteration.

## Design & Implementation
Calls `_robot_arm_control_points(q0, q_goal, pos, cfg.n_waypoints)` to prepend `q0` and append `q_goal` around the reshaped interior waypoints, then forwards to `robot_arm_hypr_path_cost_components(points, model, base_pose, obstacles, cfg; cost_cutoff)`. The keyword `cost_cutoff` defaults to `Inf`; the PSO loop passes each particle's `pbest_cost` so paths that cannot improve the personal best are reported with `total = Inf`. `pos` may be a view (`@view positions[:, pidx]`) and is not modified.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | Any | n/a | yes | Positional argument `pos`. |
| in | `cost_cutoff` | Any | n/a | no | Keyword argument `cost_cutoff` (default `Inf`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `evaluate_position`. Returns `robot_arm_hypr_path_cost_components(points, model, base_pose, obstacles, cfg; co`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:461-461`
- [[gnc.pso_path_planning_evaluate_swarm_bang|evaluate_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:409-409`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:368-368`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:268-268`

**Downstream**

- `callees` → [[gnc.config_robotarmhyprresult|RobotArmHYPRResult]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:352-352`
- `callees` → [[gnc.planner_core__robot_arm_hypr_post_refine_points|_robot_arm_hypr_post_refine_points]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:334-334`
- `callees` → [[gnc.planner_core__robot_arm_plan_from_q_reference|_robot_arm_plan_from_q_reference]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:339-339`
- `callees` → [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:270-270`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:289-289`
- `callees` → [[gnc.robot_arm_planning__reference_times|_reference_times]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:336-336`
- `callees` → [[gnc.rrt_warmstart__robot_arm_rrt_warmstart_fields|_robot_arm_rrt_warmstart_fields]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:351-351`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_control_points|_robot_arm_control_points]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:269-269`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_base_wrench_ratios|_robot_arm_hypr_base_wrench_ratios]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:340-340`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_cull_swarm_bang|_robot_arm_hypr_cull_swarm!]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:306-306`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_early_stopping_feasible|_robot_arm_hypr_early_stopping_feasible]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:291-291`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_iteration_weights|_robot_arm_hypr_iteration_weights]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:319-319`
- `callees` → [[gnc.swarm_and_retiming__robot_arm_hypr_material_improvement|_robot_arm_hypr_material_improvement]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:292-292`
- `callees` → [[gnc.swarm_and_retiming_robot_arm_sample_hypr_path|robot_arm_sample_hypr_path]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:335-335`
- `callees` → [[gncz.swarm_and_retiming__robot_arm_hypr_retime_reference|_robot_arm_hypr_retime_reference]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:338-338`
<!-- vulcan:connections:end -->

## Limitations
It captures `q0`, `q_goal`, `model`, `base_pose`, `obstacles` and `cfg` from the enclosing scope, so it cannot be tested independently. Because the cutoff in `robot_arm_hypr_path_cost_components` is applied after full evaluation, passing `pbest_cost` gives no speed benefit, only a sentinel `Inf` total. The closure is allocated every call to the planner and captures boxed variables, which can cost type stability in the hot PSO loop.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/planner_core.jl` line 268.
