---
id: gnc.path_retiming_rpo_remove_near_duplicate_samples
label: rpo_remove_near_duplicate_samples
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_retiming.jl
  symbol: rpo_remove_near_duplicate_samples
  lines:
  - 16
  - 16
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: tol
  type: Real
  units: n/a
  required: false
  description: Keyword argument `tol` (default `1.0e-10`).
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
  description: Return value of `rpo_remove_near_duplicate_samples`. Returns `pts[:,
    keep]`.
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

# rpo_remove_near_duplicate_samples

## Purpose
Filters a sampled RPO path so that no two consecutive retained samples are closer than `tol` metres, preventing zero-length segments that would otherwise produce division-by-zero or zero curvature in the retiming pipeline.

## Design & Implementation
`points` is copied into a `Matrix{Float64}`; if it has at most one column it is returned as-is. A `keep::Vector{Int}` starts with index 1, and for `j in 2:n` the column is kept only when `norm(pts[:, j] - pts[:, keep[end]]) > Float64(tol)`, so the comparison is always against the last retained sample rather than the immediately preceding raw one. When any samples are dropped, an `@warn` is emitted with `removed` and `kept` counts. The result is the column subset `pts[:, keep]`. Default `tol = 1.0e-10`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `tol` | Real | n/a | no | Keyword argument `tol` (default `1.0e-10`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_remove_near_duplicate_samples`. Returns `pts[:, keep]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:131-131`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:25-25`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:26-26`
<!-- vulcan:connections:end -->

## Limitations
The comparison against the last kept sample means a slow drift of many sub-tolerance steps can still be collapsed into a single jump larger than `tol`; this is intentional but changes path geometry. The first sample is unconditionally kept, so the final sample may be dropped if it coincides with an earlier retained point, shortening the path end. The warning fires on every call that drops samples, which can be noisy inside batch runs. Negative `tol` keeps every sample.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_retiming.jl` line 16.
