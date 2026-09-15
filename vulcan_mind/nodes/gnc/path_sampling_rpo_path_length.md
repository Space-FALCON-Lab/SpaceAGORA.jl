---
id: gnc.path_sampling_rpo_path_length
label: rpo_path_length
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_path_length
  lines:
  - 2
  - 2
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
  description: Return value of `rpo_path_length`. Returns `hypr_path_length(points)`.
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

# rpo_path_length

## Purpose
Reports the Euclidean length of an RPO waypoint polyline, the quantity sampling routines use to choose how many points a path needs.

## Design & Implementation
A one-line forward to the shared `hypr_path_length`, which sums the straight-line distances between consecutive columns of the three-by-N points matrix. Wrapping it under the `rpo_` prefix keeps the RPO planner's public surface self-contained so the HYPR core can be swapped or renamed without touching RPO callers.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_path_length`. Returns `hypr_path_length(points)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:126-126`
- [[gnc.path_sampling_rpo_sample_path_bezier|rpo_sample_path_bezier]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:25-25`
- [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:60-60`
- [[gnc.pso_refinement_rpo_refinement_shortcut_samples|rpo_refinement_shortcut_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:79-79`
- [[gnc.rrt_connect_rpo_rrt_shortcut_path|rpo_rrt_shortcut_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:236-236`
- [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:192-192`
- [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:208-208`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_path_length|hypr_path_length]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
For a Bezier control polygon this is the polygon length, not the curve's arc length, and overestimates it; callers such as `rpo_sample_path_bezier` rely on that overestimate only as a sample-count upper bound.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 2.
