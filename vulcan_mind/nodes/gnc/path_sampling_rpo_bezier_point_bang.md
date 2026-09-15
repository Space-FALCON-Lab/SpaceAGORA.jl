---
id: gnc.path_sampling_rpo_bezier_point_bang
label: rpo_bezier_point!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_bezier_point!
  lines:
  - 12
  - 12
inputs:
- id: out
  type: Any
  units: n/a
  required: true
  description: Positional argument `out`.
- id: work
  type: Any
  units: n/a
  required: true
  description: Positional argument `work`.
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `rpo_bezier_point!`; mutates `out` in place. Returns
    `hypr_bezier_point!(out, work, points, t)`.
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

# rpo_bezier_point!

## Purpose
The allocation-free Bezier evaluator used inside every sampling loop, writing the curve point into caller-owned buffers.

## Design & Implementation
Forwards to `hypr_bezier_point!` with `out` receiving the three-vector result and `work`, a matrix the same shape as `points`, serving as de Casteljau scratch space so the control polygon is never mutated. Requiring `t::Float64` avoids a conversion on each of the thousands of calls a PSO evaluation makes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `out` | Any | n/a | yes | Positional argument `out`. |
| in | `work` | Any | n/a | yes | Positional argument `work`. |
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_bezier_point!`; mutates `out` in place. Returns `hypr_bezier_point!(out, work, points, t)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_bezier_speed_estimate|rpo_bezier_speed_estimate]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:192-192`
- [[gnc.path_sampling_rpo_sample_path_bezier|rpo_sample_path_bezier]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:30-30`
- [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:213-213`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_bezier_point_bang|hypr_bezier_point!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
The caller must size `work` to match `points` exactly; a mismatch surfaces as a bounds error inside the HYPR core rather than a clear message here.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 12.
