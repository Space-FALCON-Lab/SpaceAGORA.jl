---
id: gnc.path_retiming_rpo_interpolate_along_path
label: rpo_interpolate_along_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_retiming.jl
  symbol: rpo_interpolate_along_path
  lines:
  - 39
  - 39
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
- id: s_vals
  type: Any
  units: n/a
  required: true
  description: Positional argument `s_vals`.
- id: s_query
  type: Real
  units: n/a
  required: true
  description: Positional argument `s_query`.
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
  description: Return value of `rpo_interpolate_along_path`. Returns `copy(pts[:,
    idx])` or `(1.0 - α) .* pts[:, idx] .+ α .* pts[:, idx + 1]`.
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

# rpo_interpolate_along_path

## Purpose
Returns the 3-D position on a sampled RPO path at an arbitrary arc-length coordinate `s_query`, linearly interpolating between the two bracketing samples. It is called once per propagation step inside `rpo_retime_path` to build the position reference history.

## Theory & Math
For segment $[s_i, s_{i+1}]$ with $s_i \le s_q < s_{i+1}$: $\alpha = \dfrac{s_q - s_i}{s_{i+1} - s_i}$, $\mathbf{r}(s_q) = (1-\alpha)\,\mathbf{p}_i + \alpha\,\mathbf{p}_{i+1}$, with $\alpha$ clamped to $[0,1]$.

## Design & Implementation
Inputs are coerced to `Matrix{Float64}` and `Vector{Float64}`. Three early exits clamp the query: a single sample returns column 1; `sq <= s[1]` returns the first column; `sq >= s[end]` returns the last. Otherwise `idx = clamp(searchsortedlast(s, sq), 1, n-1)` locates the segment, and a `while` loop advances `idx` past any segment whose length is `<= eps(Float64)`. If the loop reaches `n` the last column is returned. The blend factor `α = clamp((sq - s[idx]) / denom, 0, 1)` produces `(1-α)·p[idx] + α·p[idx+1]`. Every return path yields a fresh `copy` or newly allocated vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `s_vals` | Any | n/a | yes | Positional argument `s_vals`. |
| in | `s_query` | Real | n/a | yes | Positional argument `s_query`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_interpolate_along_path`. Returns `copy(pts[:, idx])` or `(1.0 - α) .* pts[:, idx] .+ α .* pts[:, idx + 1]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_resample_polyline_points|rpo_resample_polyline_points]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:52-52`
- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:250-250`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:42-42`
<!-- vulcan:connections:end -->

## Limitations
Assumes `s_vals` is monotonically non-decreasing; `searchsortedlast` gives undefined segment selection otherwise. Requires `length(s_vals) == size(points, 2)` but never checks it, so a mismatch produces a `BoundsError` or silently wrong interpolation. Interpolation is piecewise linear only, so the interpolated point can lie inside the chord of a curved path. Each call allocates full copies of `points` and `s_vals`, which is costly when invoked thousands of times per retiming.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_retiming.jl` line 39.
