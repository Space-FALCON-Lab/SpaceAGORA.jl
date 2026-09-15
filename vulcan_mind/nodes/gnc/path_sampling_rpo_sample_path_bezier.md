---
id: gnc.path_sampling_rpo_sample_path_bezier
label: rpo_sample_path_bezier
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_sample_path_bezier
  lines:
  - 17
  - 17
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
  description: Return value of `rpo_sample_path_bezier`. Returns `out`.
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

# rpo_sample_path_bezier

## Purpose
Discretises a Bezier candidate path at approximately uniform arc-length spacing so clearance and cost can be evaluated pointwise.

## Design & Implementation
Converts the control points to a dense `Matrix{Float64}` and returns immediately for a single point. It sizes the sample count as the ceiling of the larger of polygon length and end-to-end chord divided by `ds`, at least one, then evaluates the curve at `j/n` for `j` from zero to `n` into a preallocated output using the in-place evaluator and a reused work buffer. The first and last columns are overwritten with the exact control endpoints so floating-point evaluation error cannot move the path's start or goal.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `ds` | Real | n/a | yes | Positional argument `ds`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_sample_path_bezier`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:258-258`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:25-25`
- `callees` → [[gnc.path_sampling_rpo_bezier_point_bang|rpo_bezier_point!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:30-30`
- `callees` → [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:25-25`
<!-- vulcan:connections:end -->

## Limitations
Uniform steps in the parameter are not uniform in arc length — the curve moves faster where control points are far apart — so the spacing is only approximate and can exceed `ds` locally on strongly curved polygons.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 17.
