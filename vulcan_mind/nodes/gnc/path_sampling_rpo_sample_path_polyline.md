---
id: gnc.path_sampling_rpo_sample_path_polyline
label: rpo_sample_path_polyline
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_sample_path_polyline
  lines:
  - 58
  - 58
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
  description: Return value of `rpo_sample_path_polyline`. Returns `rpo_resample_polyline_points(pts,
    n)`.
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

# rpo_sample_path_polyline

## Purpose
Samples a polyline path at a requested spacing by converting the spacing into a point count and resampling.

## Design & Implementation
Measures total length with `rpo_path_length`, sets `n` to the ceiling of length over `ds` plus one with a floor of two, and delegates to `rpo_resample_polyline_points`. The plus one ensures both endpoints are included at the requested resolution.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `ds` | Real | n/a | yes | Positional argument `ds`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_sample_path_polyline`. Returns `rpo_resample_polyline_points(pts, n)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:259-259`
- [[gnc.pso_adaptive_policy_rpo_estimate_geometry_complexity|rpo_estimate_geometry_complexity]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:4-4`
- [[gnc.pso_adaptive_policy_rpo_probe_geometry_metrics|rpo_probe_geometry_metrics]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:14-14`
- [[gnc.pso_path_planning_seed_control_points|seed_control_points]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:304-304`
- [[gnc.replanning_rpo_remaining_reference_path|rpo_remaining_reference_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:193-193`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:304-304`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:61-61`
- `callees` → [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:60-60`
- `callees` → [[gnc.path_sampling_rpo_resample_polyline_points|rpo_resample_polyline_points]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:62-62`
<!-- vulcan:connections:end -->

## Limitations
Because it resamples rather than inserting points between existing vertices, the output does not pass exactly through interior waypoints; a caller needing vertex preservation must use the adaptive segment sampler instead.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 58.
