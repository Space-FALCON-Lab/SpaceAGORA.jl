---
id: gnc.path_sampling_rpo_adaptive_sampling_min_ds_m
label: rpo_adaptive_sampling_min_ds_m
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_adaptive_sampling_min_ds_m
  lines:
  - 73
  - 73
inputs:
- id: base_ds
  type: Real
  units: n/a
  required: true
  description: Positional argument `base_ds`.
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
  description: Return value of `rpo_adaptive_sampling_min_ds_m`. Returns `max(min_ds,
    1.0e-9)`.
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

# rpo_adaptive_sampling_min_ds_m

## Purpose
Computes the tightest sample spacing adaptive sampling may use near obstacles, so no path sample can straddle the keep-out volume unnoticed.

## Design & Implementation
Starts from `base_ds` floored at 1e-9 and returns it unchanged if `cfg.adaptive_sampling_enable` is off. Otherwise it reduces the minimum to the smaller of `adaptive_sampling_safe_distance_fraction` times the safety distance, when that is positive, and `adaptive_sampling_obstacle_guard_fraction` times the inflated obstacle radius, when that is positive. The final floor at 1e-9 prevents a zero step that would stall the sampling loop.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `base_ds` | Real | n/a | yes | Positional argument `base_ds`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_adaptive_sampling_min_ds_m`. Returns `max(min_ds, 1.0e-9)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_sample_path_polyline_adaptive|rpo_sample_path_polyline_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:166-166`
- [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:45-45`
- [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:206-206`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:79-79`
- `callees` → [[gnc.path_sampling_rpo_inflated_obstacle_radius_m|rpo_inflated_obstacle_radius_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:81-81`
<!-- vulcan:connections:end -->

## Limitations
Both fractions are read from configuration without range checks, so a fraction above one makes the minimum larger than intended and defeats the guard.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 73.
