---
id: gncy.pso_refinement_rpo_post_refine_path
label: rpo_post_refine_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_post_refine_path
  lines:
  - 272
  - 331
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: raw_path
  type: Matrix
  units: m
  required: true
  description: Three-row RTN path produced by a planner, to be shortened, tightened,
    and simplified without losing feasibility.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: refined_path
  type: Tuple
  units: n/a
  description: Refined path together with its final cost and a flag indicating whether
    any refinement operator was accepted.
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

# rpo_post_refine_path

## Purpose
`rpo_post_refine_path` is the shared post-processing stage that every RPO planner runs on its output. It applies three local improvement operators in rounds and keeps only changes that improve the cost, so refinement can never make a plan worse than the planner produced.

## Model & Assumptions
Refinement decisions are made under a separate configuration from the one used for the final score. `rpo_refinement_config` derives that decision configuration, which lets refinement evaluate candidates at a different sampling density or weighting than the reported cost while the returned cost is always recomputed under the caller's configuration. Acceptance is governed by `rpo_refinement_better` and `rpo_try_accept_refinement`, which compare full cost component tuples rather than scalars so a candidate that trades length for clearance can be judged on both.

## Design & Implementation
Three operators run in fixed order within each round. `rpo_refine_shortcut_refit` samples the path, finds shortcut opportunities that pass the segment safety check, and refits a Bezier through the survivors. `rpo_refine_tighten_handles` pulls control handles inward toward the path. `rpo_refine_lower_degree` attempts to represent the same path with fewer control points. Each returns a candidate, its components, and an acceptance flag; an accepted candidate immediately becomes the current path and its components become the new baseline for the operators that follow in the same round. A round that accepts nothing breaks the loop early rather than burning the remaining budget. When `refinement_enable` is false or `refinement_rounds` is zero the path is returned unchanged with its cost and a false improvement flag.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `raw_path` | Matrix | m | yes | Three-row RTN path produced by a planner, to be shortened, tightened, and simplified without losing feasibility. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `refined_path` | Tuple | n/a | — | Refined path together with its final cost and a flag indicating whether any refinement operator was accepted. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:613-613`
- [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:466-466`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:597-597`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:502-502`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:613-613`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:399-399`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:316-316`

**Downstream**

- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:275-275`
- `callees` → [[gnc.pso_refinement_rpo_refine_lower_degree|rpo_refine_lower_degree]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:312-312`
- `callees` → [[gnc.pso_refinement_rpo_refine_shortcut_refit|rpo_refine_shortcut_refit]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:284-284`
- `callees` → [[gnc.pso_refinement_rpo_refine_tighten_handles|rpo_refine_tighten_handles]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:298-298`
- `callees` → [[gnc.pso_refinement_rpo_refinement_config|rpo_refinement_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:274-274`
- `callees` → [[gncy.path_costs_rpo_path_cost|rpo_path_cost]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:277-277`
<!-- vulcan:connections:end -->

## Limitations
Because operators run in a fixed sequence and each accepts greedily, the outcome depends on operator order and the result is a local optimum rather than the best available refinement. Feasibility rests on the discrete segment sampling in `rpo_refinement_segment_is_safe`, so a thin obstacle between samples can be shortcut straight through. Degree lowering changes the curve representation, which can shift the retimed velocity profile even when the geometric cost improves. The final cost is recomputed under a different configuration than the acceptance decisions used, so an accepted candidate can score marginally worse in the returned number.

## Provenance
Mapped from pso_refinement.jl lines 272-331; include site observed at guidance_hooks.jl line 71.
