---
id: gnc.pso_refinement_rpo_refinement_segment_samples
label: rpo_refinement_segment_samples
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_segment_samples
  lines:
  - 29
  - 29
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
  description: Return value of `rpo_refinement_segment_samples`. Returns `samples`.
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

# rpo_refinement_segment_samples

## Purpose
Places evenly spaced samples along a straight segment for a collision check when adaptive sampling is disabled.

## Design & Implementation
Converts both endpoints to `SVector`, sets the segment count to the ceiling of the chord length over `ds` with `ds` floored at 1e-9 and at least one segment, and fills a three-by-`(n+1)` matrix by linear interpolation at `k/n`. Both endpoints are included exactly at `k = 0` and `k = n`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `a` | Any | n/a | yes | Positional argument `a`. |
| in | `b` | Any | n/a | yes | Positional argument `b`. |
| in | `ds` | Real | n/a | yes | Positional argument `ds`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_segment_samples`. Returns `samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:57-57`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:32-32`
<!-- vulcan:connections:end -->

## Limitations
Uniform spacing is blind to where the station is, so with a coarse `ds` a segment can pass close to the keep-out between samples and still be judged safe.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 29.
