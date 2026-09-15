---
id: gnc.path_sampling_rpo_sample_path
label: rpo_sample_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_sample_path
  lines:
  - 257
  - 257
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: ds
  type: Real
  units: n/a
  required: true
  description: Positional argument `ds`.
- id: curve_type
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `curve_type` (default `:bezier`).
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
  description: Return value of `rpo_sample_path`.
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

# rpo_sample_path

## Purpose
The single entry point for sampling an RPO candidate path, choosing between polyline and Bezier and between fixed and adaptive spacing from configuration.

## Design & Implementation
Two methods. The `(points, ds; curve_type)` form dispatches to the fixed-spacing Bezier or polyline sampler. The `(points, cfg, geometry; ...)` form reads `safe_distance_m`, `sample_ds_m` and `curve_type` from the `RPOPSOConfig` and, when `adaptive_sampling_enable` is set, routes to the adaptive Bezier or polyline sampler with the station geometry; otherwise it falls back to the fixed form. Any other `curve_type` raises `ArgumentError` naming the two supported symbols.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `ds` | Real | n/a | yes | Positional argument `ds`. |
| in | `curve_type` | Symbol | n/a | no | Keyword argument `curve_type` (default `:bezier`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_sample_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:93-93`
- [[gnc.pso_path_planning_seed_control_points|seed_control_points]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:303-303`
- [[gnc.pso_refinement_rpo_refine_lower_degree|rpo_refine_lower_degree]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:236-236`
- [[gnc.pso_refinement_rpo_refine_shortcut_refit|rpo_refine_shortcut_refit]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:175-175`
- [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:442-442`
- [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:183-183`
- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:123-123`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:303-303`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_sample_path_bezier|rpo_sample_path_bezier]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:258-258`
- `callees` → [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:259-259`
- `callees` → [[gnc.path_sampling_rpo_sample_path_polyline_adaptive|rpo_sample_path_polyline_adaptive]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:282-282`
- `callees` → [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:275-275`
<!-- vulcan:connections:end -->

## Limitations
The fixed-spacing form has no geometry argument, so a caller who forgets to pass `geometry` to the configured form gets fixed sampling even with adaptive enabled in `cfg` only if they call the wrong method — the two signatures are easy to confuse because both accept `points` first.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 257.
