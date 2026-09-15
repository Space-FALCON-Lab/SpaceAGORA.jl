---
id: gnc.hypr_utils_hypr_bezier_point
label: hypr_bezier_point
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_bezier_point
  lines:
  - 28
  - 28
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
  description: Return value of `hypr_bezier_point`. Returns `hypr_bezier_point!(out,
    work, pts, Float64(t))`.
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

# hypr_bezier_point

## Purpose
Allocating convenience wrapper that evaluates a Bezier control polygon at one parameter value.

## Design & Implementation
Converts the points, allocates an output vector of the row dimension and a work matrix the same shape as the points, converts `t` to `Float64`, and delegates to the in-place evaluator. Because it allocates its own buffers, it never disturbs caller-owned work space.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_bezier_point`. Returns `hypr_bezier_point!(out, work, pts, Float64(t))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_bezier_point|rpo_bezier_point]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:8-8`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:32-32`
- `callees` → [[gnc.hypr_utils_hypr_bezier_point_bang|hypr_bezier_point!]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:32-32`
<!-- vulcan:connections:end -->

## Limitations
Two allocations per call, so it must not be used inside sampling loops; that is what the `!` variant exists for.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 28.
