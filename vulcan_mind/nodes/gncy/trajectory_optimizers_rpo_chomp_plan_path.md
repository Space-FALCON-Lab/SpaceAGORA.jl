---
id: gncy.trajectory_optimizers_rpo_chomp_plan_path
label: rpo_chomp_plan_path
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_chomp_plan_path
  lines:
  - 225
  - 333
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: chomp_request
  type: Tuple
  units: n/a
  required: true
  description: Start and goal RTN positions, station geometry, and the PSO configuration
    whose weights and waypoint count seed the optimizer.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: chomp_result
  type: NamedTuple
  units: n/a
  description: Refined path, raw path, refined and raw costs, cost components, iteration
    count, and cost history for the CHOMP comparison run.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# rpo_chomp_plan_path

## Purpose
`rpo_chomp_plan_path` implements the CHOMP-style gradient comparison planner used as a baseline against the HYPR PSO planner. It optimises the interior waypoints of a start-to-goal trajectory under a smoothness-plus-obstacle objective using a covariant gradient step.

## Theory & Math
The covariant update is $\theta \leftarrow \theta - \eta\, R^{-1} \nabla_\theta J$, where $R$ is the second-difference smoothness metric and $J$ combines the smoothness cost with the soft obstacle potential. Central differences give $\partial J/\partial \theta_{d,i} \approx (J(\theta + h e_{d,i}) - J(\theta - h e_{d,i}))/(2h)$.

## Model & Assumptions
The trajectory is parameterised only by its interior waypoints, with endpoints pinned, so the optimiser cannot drift the start or the goal. The objective combines a soft obstacle potential with a second-difference smoothness term weighted by `settings.w_smooth`. The obstacle weight is scaled by `optimizer.w_obs_scale` relative to the shared PSO configuration so the comparison planner can be tuned without changing the reference cost used for scoring. The effective obstacle margin is the larger of the optimizer margin and the safe distance, which guarantees the comparison never plans inside the keep-out shell the scorer enforces.

## Design & Implementation
Interior waypoints are seeded from `initial_path` when a warm start is supplied and otherwise interpolated. Search bounds come from `rpo_trajectory_search_bounds` and are then relaxed to contain the seed so a warm start is never clipped on the first step. The covariant metric inverse `R_inv` from `rpo_second_difference_metric` premultiplies each spatial row of the numeric gradient, which is itself formed by central differences with step size `max(gradient_eps, 1e-6 * max(abs(theta), 1))`. Each iteration performs a backtracking line search over the scale sequence 1.0, 0.5, 0.25, 0.125, 0.0625 and accepts the first strictly improving trial. Termination is triggered by runtime limit, a rejected step, or `no_change_iters` consecutive iterations within `no_change_tol`. The accepted path is finally passed through `rpo_post_refine_path` so the comparison shares the same refinement stage as HYPR.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `chomp_request` | Tuple | n/a | yes | Start and goal RTN positions, station geometry, and the PSO configuration whose weights and waypoint count seed the optimizer. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `chomp_result` | NamedTuple | n/a | — | Refined path, raw path, refined and raw costs, cost components, iteration count, and cost history for the CHOMP comparison run. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:380-380`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:245-245`
- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:258-258`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:302-302`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:236-236`
- `callees` → [[gnc.trajectory_optimizers_rpo_chomp_numeric_gradient|rpo_chomp_numeric_gradient]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:273-273`
- `callees` → [[gnc.trajectory_optimizers_rpo_clamp_internal_waypoints|rpo_clamp_internal_waypoints]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:283-283`
- `callees` → [[gnc.trajectory_optimizers_rpo_second_difference_metric|rpo_second_difference_metric]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:244-244`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_internal_points_from_seed|rpo_trajectory_internal_points_from_seed]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:240-240`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_points_from_internal|rpo_trajectory_points_from_internal]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:256-256`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_search_bounds|rpo_trajectory_search_bounds]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:241-241`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:247-247`
- `callees` → [[gnc.trajectory_optimizers_rpochompsettings|RPOCHOMPSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:231-231`
- `callees` → [[gnc.trajectory_optimizers_rpotrajectoryoptimizersettings|RPOTrajectoryOptimizerSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:232-232`
- `callees` → [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:316-316`
<!-- vulcan:connections:end -->

## Limitations
The gradient is numeric, so each iteration costs six objective evaluations per interior waypoint and scales poorly with waypoint count. The soft obstacle potential admits paths that violate the hard keep-out constraint; feasibility is only established later by the scorer and the tracker. The line search accepts the first improving scale rather than the best, and if all five scales fail the iteration stops even when a smaller step would have helped.

## Provenance
Mapped from trajectory_optimizers.jl lines 225-333; include site observed at guidance_hooks.jl line 75.
