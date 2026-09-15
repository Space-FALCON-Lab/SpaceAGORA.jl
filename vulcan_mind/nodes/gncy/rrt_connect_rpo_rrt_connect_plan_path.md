---
id: gncy.rrt_connect_rpo_rrt_connect_plan_path
label: rpo_rrt_connect_plan_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_connect_plan_path
  lines:
  - 308
  - 418
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: rrt_request
  type: Tuple
  units: m
  required: true
  description: Start and goal RTN positions, station geometry, PSO configuration,
    and RRT-Connect settings governing iteration count and step size.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: rrt_result
  type: NamedTuple
  units: n/a
  description: Refined and raw paths with their costs and components, iteration count,
    and a flag stating whether a connected path was actually found.
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

# rpo_rrt_connect_plan_path

## Purpose
`rpo_rrt_connect_plan_path` is the bidirectional sampling-based planner used both as a comparison baseline and as the warm-start generator for the PSO planner. It grows a tree from the start and a tree from the goal, alternating which tree extends toward a random sample and which attempts to connect to it.

## Model & Assumptions
The planner works on polyline paths, so the configuration is locally rebuilt with `curve_type = :polyline` before any cost is evaluated; this keeps the tree geometry and the scoring geometry consistent. The straight-line case is tested first, and if `rpo_rrt_segment_is_safe` clears the direct segment the function returns immediately with the two-point path, avoiding tree construction entirely for unobstructed transfers. Sampling is confined to the same bounds the PSO planner uses, obtained from `rpo_pso_bounds`, so warm starts are guaranteed to lie inside the swarm search box.

## Design & Implementation
Each iteration picks the growing tree by parity of the iteration index, so the roles alternate strictly. The random sample is replaced by the opposite tree root with probability `goal_sample_rate`, clamped to the unit interval, which biases growth toward connection. `rpo_rrt_extend!` adds one step toward the sample and returns a status; a `:trapped` status skips to the next iteration. Otherwise `rpo_rrt_connect!` drives the other tree repeatedly toward the newly added node, and a `:reached` status joins both branches through `rpo_rrt_connect_join_paths`. A runtime budget is checked from the second iteration onward. The raw path is then shortened by `rpo_rrt_shortcut_path` and passed through `rpo_post_refine_path`, and both raw and refined costs are reported so the value of post-processing is measurable.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `rrt_request` | Tuple | m | yes | Start and goal RTN positions, station geometry, PSO configuration, and RRT-Connect settings governing iteration count and step size. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `rrt_result` | NamedTuple | n/a | — | Refined and raw paths with their costs and components, iteration count, and a flag stating whether a connected path was actually found. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_rrt_connect_seed_path|_rpo_comparison_rrt_connect_seed_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:244-244`
- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:302-302`
- [[gnc.pso_path_planning_rpo_pso_rrt_warmstart_path|rpo_pso_rrt_warmstart_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:134-134`
- [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:431-431`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:325-325`
- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:334-334`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:318-318`
- `callees` → [[gnc.rrt_connect_rpo_rrt_connect_bang|rpo_rrt_connect!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:380-380`
- `callees` → [[gnc.rrt_connect_rpo_rrt_connect_join_paths|rpo_rrt_connect_join_paths]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:390-390`
- `callees` → [[gnc.rrt_connect_rpo_rrt_extend_bang|rpo_rrt_extend!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:370-370`
- `callees` → [[gnc.rrt_connect_rpo_rrt_random_state|rpo_rrt_random_state]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:366-366`
- `callees` → [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:327-327`
- `callees` → [[gnc.rrt_connect_rpo_rrt_shortcut_path|rpo_rrt_shortcut_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:398-398`
- `callees` → [[gnc.rrt_connect_rporrtconnectsettings|RPORRTConnectSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:314-314`
- `callees` → [[gnc.rrt_connect_rporrtconnecttree|RPORRTConnectTree]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:353-353`
- `callees` → [[gncy.pso_helpers_rpo_pso_bounds|rpo_pso_bounds]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:321-321`
- `callees` → [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:399-399`
<!-- vulcan:connections:end -->

## Limitations
When no connection is found the function returns the direct start-to-goal segment with `path_found` set false, so a caller that ignores that flag will silently fly a path already known to be unsafe. Sampling is uniform over the axis-aligned box with no obstacle awareness, so narrow passages need many iterations. The step size and collision resolution come from the settings record, and a step larger than the collision sampling resolution can tunnel through thin geometry. Reported iteration count reflects extend attempts rather than nodes added.

## Provenance
Mapped from rrt_connect.jl lines 308-418; include site observed at guidance_hooks.jl line 72.
