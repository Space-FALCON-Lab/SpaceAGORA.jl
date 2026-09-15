---
id: gnc.path_sampling_rpo_resample_polyline_points
label: rpo_resample_polyline_points
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_resample_polyline_points
  lines:
  - 45
  - 45
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: n_samples
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_samples`.
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
  description: Return value of `rpo_resample_polyline_points`. Returns `out`.
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

# rpo_resample_polyline_points

## Purpose
Redistributes a polyline onto a fixed number of points equally spaced by arc length, used when downstream code needs a known column count.

## Design & Implementation
Forces at least two samples, computes cumulative arc-length parameters with `rpo_arc_length_params`, then places sample `j` at fraction `(j-1)/(n_samples-1)` of the total length through `rpo_interpolate_along_path`. Because the spacing is computed from the cumulative parameter rather than the segment index, the original vertices need not be evenly spaced.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `n_samples` | Int | n/a | yes | Positional argument `n_samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_resample_polyline_points`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:62-62`
- [[gnc.pso_path_planning_seed_control_points|seed_control_points]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:307-307`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:307-307`

**Downstream**

- `callees` → [[gnc.path_retiming_rpo_arc_length_params|rpo_arc_length_params]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:48-48`
- `callees` → [[gnc.path_retiming_rpo_interpolate_along_path|rpo_interpolate_along_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:52-52`
<!-- vulcan:connections:end -->

## Limitations
Original vertices are generally not preserved except the two endpoints, so a sharp corner in the input can be rounded off by interpolation between samples that straddle it.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 45.
