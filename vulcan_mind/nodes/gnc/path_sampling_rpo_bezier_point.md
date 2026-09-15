---
id: gnc.path_sampling_rpo_bezier_point
label: rpo_bezier_point
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_bezier_point
  lines:
  - 7
  - 7
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: t
  type: Real
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
  description: Return value of `rpo_bezier_point`. Returns `hypr_bezier_point(points,
    t)`.
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

# rpo_bezier_point

## Purpose
Evaluates the Bezier curve defined by a control polygon at a normalised parameter, the allocating convenience form.

## Design & Implementation
Forwards `points` and `t` to `hypr_bezier_point`, which applies de Casteljau subdivision over the columns of the control matrix and returns a fresh three-vector. Being allocation-per-call, it is intended for one-off evaluations such as plotting or endpoint checks rather than inner sampling loops, which use the in-place variant.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_bezier_point`. Returns `hypr_bezier_point(points, t)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_bezier_point|hypr_bezier_point]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:8-8`
<!-- vulcan:connections:end -->

## Limitations
`t` is not clamped to the unit interval, so an out-of-range parameter extrapolates the polynomial beyond the control polygon without error.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 7.
