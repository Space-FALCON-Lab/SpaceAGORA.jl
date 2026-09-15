---
id: gnc.path_sampling_rpo_adaptive_segment_samples
label: rpo_adaptive_segment_samples
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_adaptive_segment_samples
  lines:
  - 110
  - 110
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
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
- id: min_ds_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `min_ds_m`.
- id: max_ds_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `max_ds_m`.
- id: far_clearance_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `far_clearance_m`.
- id: power
  type: Real
  units: n/a
  required: false
  description: Keyword argument `power` (default `1.0`).
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
  description: Return value of `rpo_adaptive_segment_samples`. Returns `out`.
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

# rpo_adaptive_segment_samples

## Purpose
Samples one straight segment between two waypoints with clearance-adaptive spacing while guaranteeing both endpoints appear exactly.

## Design & Implementation
Returns the single start point for a degenerate segment shorter than machine epsilon. Otherwise it walks from `a` toward `b` along the unit direction, at each position querying `rpo_clearance_distance_to_station` and asking `rpo_adaptive_sampling_step_m` for the next step, appending points until the distance is covered. A `max_steps` bound of `dist / min_ds + 2` guards against a runaway loop, and the final sample is overwritten with `b` exactly. Points are gathered as `SVector`s and copied into a dense matrix at the end.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Any | n/a | yes | Positional argument `a`. |
| in | `b` | Any | n/a | yes | Positional argument `b`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `min_ds_m` | Real | n/a | yes | Keyword argument `min_ds_m`. |
| in | `max_ds_m` | Real | n/a | yes | Keyword argument `max_ds_m`. |
| in | `far_clearance_m` | Real | n/a | yes | Keyword argument `far_clearance_m`. |
| in | `power` | Real | n/a | no | Keyword argument `power` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_adaptive_segment_samples`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_sample_path_polyline_adaptive|rpo_sample_path_polyline_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:170-170`
- [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:46-46`
- [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:109-109`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:130-130`
- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:134-134`
- `callees` → [[gnc.path_sampling_rpo_adaptive_sampling_step_m|rpo_adaptive_sampling_step_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:135-135`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:128-128`
<!-- vulcan:connections:end -->

## Limitations
Clearance is evaluated at the current sample, not along the step, so a segment that grazes the station between two samples can still be under-resolved if the step cap is loose; hitting `max_steps` silently truncates rather than warning.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 110.
