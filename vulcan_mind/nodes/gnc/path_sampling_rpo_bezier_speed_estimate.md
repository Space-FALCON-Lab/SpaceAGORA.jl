---
id: gnc.path_sampling_rpo_bezier_speed_estimate
label: rpo_bezier_speed_estimate
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_bezier_speed_estimate
  lines:
  - 187
  - 187
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: work
  type: Any
  units: n/a
  required: true
  description: Positional argument `work`.
- id: point
  type: Any
  units: n/a
  required: true
  description: Positional argument `point`.
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
  description: Return value of `rpo_bezier_speed_estimate`. Returns `norm(p2 - point)
    / abs(t2 - t)`.
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

# rpo_bezier_speed_estimate

## Purpose
Estimates the magnitude of the Bezier curve's parametric derivative at a point, so adaptive sampling can convert a desired spatial step into a parameter increment.

## Design & Implementation
Chooses a probe parameter 1e-4 ahead of `t`, or behind it if already at the end, and returns zero if neither moves. It evaluates the curve at the probe into a temporary through the in-place evaluator and returns the chord length divided by the parameter difference — a forward finite difference of the position with respect to `t`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `work` | Any | n/a | yes | Positional argument `work`. |
| in | `point` | Any | n/a | yes | Positional argument `point`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_bezier_speed_estimate`. Returns `norm(p2 - point) / abs(t2 - t)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:227-227`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_bezier_point_bang|rpo_bezier_point!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:192-192`
<!-- vulcan:connections:end -->

## Limitations
A fixed 1e-4 probe step is not scaled to the curve's length, so on very short curves the finite difference is dominated by rounding while on very long ones it is coarse; the temporary `p2` is allocated on every call.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 187.
