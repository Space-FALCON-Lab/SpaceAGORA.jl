---
id: gncy.pso_path_planning_rpo_pso_plan_path
label: rpo_pso_plan_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_plan_path
  lines:
  - 184
  - 631
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: plan_request
  type: Tuple
  units: m
  required: true
  description: Start and goal RTN positions, station geometry, and the base PSO configuration
    for the transfer to be planned.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: plan_result
  type: NamedTuple
  units: n/a
  description: Refined path, cost, cost components, resolved configuration, adaptive
    diagnostics, warm-start diagnostics, early-stop and timeout status, and the cost
    history.
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

# rpo_pso_plan_path

## Purpose
`rpo_pso_plan_path` is the primary RPO path planner. It sizes the search to the problem, optionally warm-starts from an RRT-Connect solution, runs a bounded particle swarm over the interior waypoints, and hands the swarm best to the refinement stage before returning a fully diagnosed plan.

## Theory & Math
The swarm update is $v_{d,i} \leftarrow w v_{d,i} + c_1 r_1 (p_{d,i} - x_{d,i}) + c_2 r_2 (g_d - x_{d,i})$ followed by $x_{d,i} \leftarrow \mathrm{clamp}(x_{d,i} + v_{d,i}, lo_d, hi_d)$, with $r_1, r_2$ drawn independently per dimension from the thread-local generator.

## Model & Assumptions
Safe distance resolution happens first and deliberately treats an explicit zero as an absent value when the configuration carries a positive safe distance, so a caller cannot accidentally disable the keep-out shell by passing a default. The configuration is then re-derived twice, once for the adaptive policy and once for the effective safe distance, so the configuration returned in the result is exactly the one the swarm ran under. The degenerate zero-waypoint case is handled analytically as a straight line rather than by running a swarm over an empty search space.

## Design & Implementation
When waypoints exist, `rpo_pso_rrt_warmstart_path` may produce a seed path. For Bezier curves the waypoint count is raised to match the warm-start path, bounded above by `reexplore_max_waypoints` when that is positive, and the configuration is rebuilt at the new count. The swarm state is preallocated as dense matrices for positions, velocities, personal bests, personal best costs and obstacle scores, plus per-particle stagnation counters. The update is the standard PSO recurrence with per-iteration weights from `rpo_pso_iteration_weights` and per-thread random generators indexed by `threadid()`, with positions clamped to the replicated bounds. Every phase checks an iteration runtime budget and records a timeout event naming the phase, either `:velocity_update` or `:swarm_evaluate`, before breaking. Finally the global best is expanded to a full path by `rpo_position_to_path`, refined by `rpo_post_refine_path`, and rescored.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `plan_request` | Tuple | m | yes | Start and goal RTN positions, station geometry, and the base PSO configuration for the transfer to be planned. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `plan_result` | NamedTuple | n/a | — | Refined path, cost, cost components, resolved configuration, adaptive diagnostics, warm-start diagnostics, early-stop and timeout status, and the cost history. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:276-276`
- [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:23-23`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:186-186`
- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:210-210`
- `callees` → [[gnc.path_sampling_rpo_resample_polyline_points|rpo_resample_polyline_points]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:307-307`
- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:303-303`
- `callees` → [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:304-304`
- `callees` → [[gnc.planner_core_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:368-368`
- `callees` → [[gnc.pso_helpers_rpo_position_to_path|rpo_position_to_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:369-369`
- `callees` → [[gnc.pso_helpers_rpo_pso_warmstart_bounds|rpo_pso_warmstart_bounds]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:326-326`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:266-266`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:189-189`
- `callees` → [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:440-440`
- `callees` → [[gnc.pso_path_planning_cfg_for_current|cfg_for_current]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:285-285`
- `callees` → [[gnc.pso_path_planning_evaluate_position|evaluate_position]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:368-368`
- `callees` → [[gnc.pso_path_planning_evaluate_swarm_bang|evaluate_swarm!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:402-402`
- `callees` → [[gnc.pso_path_planning_iteration_elapsed_s|iteration_elapsed_s]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:254-254`
- `callees` → [[gnc.pso_path_planning_iteration_timed_out|iteration_timed_out]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:259-259`
- `callees` → [[gnc.pso_path_planning_record_iteration_bang|record_iteration!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:265-265`
- `callees` → [[gnc.pso_path_planning_record_iteration_timeout_bang|record_iteration_timeout!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:272-272`
- `callees` → [[gnc.pso_path_planning_reexplore_waypoint_count|reexplore_waypoint_count]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:431-431`
- `callees` → [[gnc.pso_path_planning_reset_swarm_bang|reset_swarm!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:315-315`
- `callees` → [[gnc.pso_path_planning_rpo_pso_cull_swarm_bang|rpo_pso_cull_swarm!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:546-546`
- `callees` → [[gnc.pso_path_planning_rpo_pso_early_stopping_feasible|rpo_pso_early_stopping_feasible]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:575-575`
- `callees` → [[gnc.pso_path_planning_rpo_pso_effective_safe_distance|rpo_pso_effective_safe_distance]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:191-191`
- `callees` → [[gnc.pso_path_planning_rpo_pso_empty_warmstart_diagnostics|rpo_pso_empty_warmstart_diagnostics]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:199-199`
- `callees` → [[gnc.pso_path_planning_rpo_pso_iteration_weights|rpo_pso_iteration_weights]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:588-588`
- `callees` → [[gnc.pso_path_planning_rpo_pso_material_improvement|rpo_pso_material_improvement]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:569-569`
- `callees` → [[gnc.pso_path_planning_rpo_pso_protected_particle_mask|rpo_pso_protected_particle_mask]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:444-444`
- `callees` → [[gnc.pso_path_planning_rpo_pso_rrt_warmstart_path|rpo_pso_rrt_warmstart_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:198-198`
- `callees` → [[gnc.pso_path_planning_rpo_pso_stagnation_count_after_learning|rpo_pso_stagnation_count_after_learning]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:493-493`
- `callees` → [[gnc.pso_path_planning_seed_control_points|seed_control_points]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:290-290`
- `callees` → [[gnc.pso_path_planning_update_swarm_best_from_components_bang|update_swarm_best_from_components!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:380-380`
- `callees` → [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:305-305`
- `callees` → [[gncy.pso_adaptive_policy_rpo_adaptive_pso_config|rpo_adaptive_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:190-190`
- `callees` → [[gncy.pso_helpers_rpo_pso_bounds|rpo_pso_bounds]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:327-327`
- `callees` → [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:613-613`
<!-- vulcan:connections:end -->

## Limitations
Per-thread generators are indexed by thread id, so a run is only reproducible for a fixed thread count and scheduling pattern even with a fixed seed. The runtime budget terminates mid-iteration, which means a timed-out plan can return a global best that was never re-evaluated under the final configuration. Bounds are static for the whole run apart from the warm-start rebuild, so a swarm that converges to a box face cannot expand past it. The obstacle penalty remains soft throughout, so the returned path is not guaranteed feasible without checking the reported minimum clearance.

## Provenance
Mapped from pso_path_planning.jl lines 184-631; include site observed at guidance_hooks.jl line 73.
