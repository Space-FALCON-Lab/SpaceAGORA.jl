---
id: simulation.vacuum_predicted_gram__eval_natural_cubic_spline
label: _eval_natural_cubic_spline
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _eval_natural_cubic_spline
  lines:
  - 138
  - 138
inputs:
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: t0
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t0`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
- id: ys
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `ys`.
- id: Ms
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `Ms`.
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
  type: Float64
  units: n/a
  description: Return value of `_eval_natural_cubic_spline`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _eval_natural_cubic_spline

## Purpose
Evaluates a natural cubic spline at time `t` given uniformly spaced knot values and the second-derivative coefficients produced by `_natural_cubic_spline_build!`. `_query_vacuum_gram_cache!` uses it to reconstruct `log(ρ)` and temperature between cache knots.

## Theory & Math
On segment $i$ with $\Delta = t - t_i$, the spline is
$$S(t) = y_i + b_i\Delta + c_i\Delta^2 + d_i\Delta^3,\quad b_i = \frac{y_{i+1} - y_i}{h} - \frac{h(M_{i+1} + 2M_i)}{6},\; c_i = \frac{M_i}{2},\; d_i = \frac{M_{i+1} - M_i}{6h},$$
where $M_i$ are the natural-spline second derivatives at the knots.

## Design & Implementation
Signature `_eval_natural_cubic_spline(t, t0, h::Float64, ys::Vector{Float64}, Ms::Vector{Float64})::Float64`, `@inline`. The segment index is `idx = clamp(floor(Int, (t - t0)/h) + 1, 1, n - 1)` with `n = length(ys)`, and the local offset is `dx = t - (t0 + (idx-1) h)`. It loads `ys[idx]`, `ys[idx+1]`, `Ms[idx]`, `Ms[idx+1]` with `@inbounds`, forms the standard polynomial coefficients `b = (y_{i+1} - y_i)/h - h (M_{i+1} + 2 M_i)/6`, `c = M_i/2`, `d = (M_{i+1} - M_i)/(6h)`, and returns `y_i + dx (b + dx (c + dx d))` by Horner evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `t0` | Float64 | n/a | yes | Positional argument `t0`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `ys` | Vector{Float64} | n/a | yes | Positional argument `ys`. |
| in | `Ms` | Vector{Float64} | n/a | yes | Positional argument `Ms`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_eval_natural_cubic_spline`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:253-253`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `idx` is clamped, times outside `[t0, t0 + (n-1)h]` are extrapolated with the end cubic rather than rejected; the caller must enforce the horizon check. `@inbounds` makes a mismatch between `length(ys)` and `length(Ms)` undefined behaviour instead of an error. With `n < 2` the clamp yields `idx = 1` and `ys[2]` is read out of bounds. `h <= 0` gives division by zero.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 138.
