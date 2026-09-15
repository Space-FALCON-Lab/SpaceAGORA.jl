---
id: gnc.pso_refinement_rpo_refinement_segment_is_safe
label: rpo_refinement_segment_is_safe
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_segment_is_safe
  lines:
  - 42
  - 42
inputs:
- id: a
  type: Any
  units: n/a
  required: true
  description: Positional argument `a`.
- id: b
  type: Any
  units: n/a
  required: true
  description: Positional argument `b`.
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
  description: Return value of `rpo_refinement_segment_is_safe`. Returns `stats.min_clearance
    + 1.0e-9 >= required_clearance`.
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

# rpo_refinement_segment_is_safe

## Purpose
Decides whether a proposed straight shortcut between two path points stays clear of the station with the required margin.

## Design & Implementation
Chooses the sampling density with `rpo_hypr_sampling_density_m`, then samples the segment either adaptively — computing the minimum spacing from the geometry and using `rpo_adaptive_segment_samples` — or uniformly. The required clearance is `safe_distance_m` plus `refinement_straight_clearance_margin_m`, and the segment is accepted if the minimum clearance from `rpo_clearance_stats_from_samples` meets it within 1e-9.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Any | n/a | yes | Positional argument `a`. |
| in | `b` | Any | n/a | yes | Positional argument `b`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_segment_is_safe`. Returns `stats.min_clearance + 1.0e-9 >= required_clearance`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_refinement_shortcut_samples|rpo_refinement_shortcut_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:82-82`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:59-59`
- `callees` → [[gnc.path_costs_rpo_clearance_stats_from_samples|rpo_clearance_stats_from_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:60-60`
- `callees` → [[gnc.path_sampling_rpo_adaptive_sampling_min_ds_m|rpo_adaptive_sampling_min_ds_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:45-45`
- `callees` → [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:46-46`
- `callees` → [[gnc.pso_parameters_rpo_hypr_sampling_density_m|rpo_hypr_sampling_density_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:43-43`
- `callees` → [[gnc.pso_refinement_rpo_refinement_segment_samples|rpo_refinement_segment_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:57-57`
<!-- vulcan:connections:end -->

## Limitations
The extra straight-segment margin is applied only here, so a shortcut can be rejected even though the same geometry would pass the ordinary path cost; this is intentional conservatism but makes the two checks disagree at the boundary.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 42.
