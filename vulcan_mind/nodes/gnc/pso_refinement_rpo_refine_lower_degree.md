---
id: gnc.pso_refinement_rpo_refine_lower_degree
label: rpo_refine_lower_degree
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refine_lower_degree
  lines:
  - 232
  - 232
inputs:
- id: path
  type: Any
  units: n/a
  required: true
  description: Positional argument `path`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: current_components
  type: Any
  units: n/a
  required: true
  description: Positional argument `current_components`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
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
  description: Return value of `rpo_refine_lower_degree`. Returns `current, current_components,
    improved`.
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

# rpo_refine_lower_degree

## Purpose
One refinement strategy: try to represent the current path with fewer Bezier control points, accepting any lower degree that does not worsen the objective.

## Design & Implementation
Returns unchanged for three or fewer columns. It samples the current path densely once, then for every control count from one less than the current down to two, fits a polygon to the dense samples and submits it through `rpo_try_accept_refinement`. On acceptance it adopts the smaller polygon and resamples so later, even smaller fits are made against the new curve.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | Any | n/a | yes | Positional argument `path`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `current_components` | Any | n/a | yes | Positional argument `current_components`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refine_lower_degree`. Returns `current, current_components, improved`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:312-312`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:236-236`
- `callees` → [[gnc.pso_parameters_rpo_hypr_refinement_sampling_density_m|rpo_hypr_refinement_sampling_density_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:241-241`
- `callees` → [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:246-246`
- `callees` → [[gnc.pso_refinement_rpo_try_accept_refinement|rpo_try_accept_refinement]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:247-247`
<!-- vulcan:connections:end -->

## Limitations
The loop keeps descending after a rejection, so it can accept a much lower degree after skipping several intermediate ones, which occasionally produces a coarser path than a greedy stop-at-first-failure would; each candidate costs a fit plus a full path evaluation.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 232.
