---
id: gnc.path_sampling_rpo_sample_path_polyline_adaptive
label: rpo_sample_path_polyline_adaptive
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_sample_path_polyline_adaptive
  lines:
  - 157
  - 157
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
- id: base_ds_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `base_ds_m` (default `cfg.sample_ds_m`).
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
  description: Return value of `rpo_sample_path_polyline_adaptive`. Returns `hcat(samples...)`.
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

# rpo_sample_path_polyline_adaptive

## Purpose
Applies adaptive segment sampling across an entire polyline, producing one continuous point matrix with no duplicated interior vertices.

## Design & Implementation
Computes the minimum spacing through `rpo_adaptive_sampling_min_ds_m` and takes the maximum from `cfg.adaptive_sampling_max_ds_m` floored at that minimum. Each consecutive vertex pair is sampled by `rpo_adaptive_segment_samples` with the configured far clearance and power; for every segment after the first, the leading column is dropped because it duplicates the previous segment's final column. The segments are joined with `hcat`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `base_ds_m` | Real | n/a | no | Keyword argument `base_ds_m` (default `cfg.sample_ds_m`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_sample_path_polyline_adaptive`. Returns `hcat(samples...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:282-282`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_adaptive_sampling_min_ds_m|rpo_adaptive_sampling_min_ds_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:166-166`
- `callees` → [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:170-170`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:181-181`
<!-- vulcan:connections:end -->

## Limitations
Spacing is decided per segment independently, so the step size can jump discontinuously at a vertex; `hcat` over a splatted vector of matrices allocates once per segment plus the final concatenation.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 157.
