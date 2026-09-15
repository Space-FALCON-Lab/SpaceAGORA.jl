---
id: gnc.hypr_utils_hypr_bezier_point_bang
label: hypr_bezier_point!
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_bezier_point!
  lines:
  - 36
  - 36
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
  description: Return value of `hypr_bezier_point!`; mutates `out` in place. Returns
    `out`.
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

# hypr_bezier_point!

## Purpose
Evaluates a Bezier curve by de Casteljau's algorithm into caller-provided buffers, the allocation-free core all path sampling runs on.

## Theory & Math
$$
P^{(r)}_j = (1 - t)\, P^{(r-1)}_j + t\, P^{(r-1)}_{j+1},\qquad r = 1 \ldots n-1,\quad B(t) = P^{(n-1)}_1
$$

## Design & Implementation
Copies the `n` control columns into `work`, then for each of `n - 1` rounds replaces column `j` with the linear blend `(1-t) work[j] + t work[j+1]` for `j` up to `n - r`, under `@inbounds`. After the last round column one holds the curve point, which is copied into `out` and returned. De Casteljau was chosen over the Bernstein sum because it is numerically stable for any `t` in the unit interval and needs no binomial coefficients.

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
| out | `result` | Any | n/a | — | Return value of `hypr_bezier_point!`; mutates `out` in place. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.hypr_utils_hypr_bezier_point|hypr_bezier_point]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:32-32`
- [[gnc.path_sampling_rpo_bezier_point_bang|rpo_bezier_point!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:13-13`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`work` must have at least `n` columns and the same row count as `points`; no check is made. The broadcast on column views still allocates a temporary per blend in some Julia versions, so it is allocation-light rather than allocation-free.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 36.
