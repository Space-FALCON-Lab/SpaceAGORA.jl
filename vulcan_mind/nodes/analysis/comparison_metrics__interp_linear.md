---
id: analysis.comparison_metrics__interp_linear
label: _interp_linear
kind: function
source:
  file: src/analysis/verification/telemetry_verification/comparison_metrics.jl
  symbol: _interp_linear
  lines:
  - 1
  - 1
inputs:
- id: x
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `x`.
- id: y
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `y`.
- id: xq
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `xq`.
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
  description: Return value of `_interp_linear`. Returns `out`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _interp_linear

## Purpose

Resamples a simulated series onto the telemetry sample points so that simulation and flight data can be differenced element by element. Every comparison in the telemetry verification pipeline passes through this piecewise-linear interpolator before any error metric is computed.

## Design & Implementation

`_interp_linear(x, y, xq)` validates that `x` and `y` have equal length and that the domain is non-empty, throwing `ArgumentError` otherwise, and short-circuits a single-sample domain by returning `fill(y[1], length(xq))`. It preallocates `out` and walks the query points with a persistent bracket cursor `j` that only advances forward, giving linear rather than logarithmic cost per point when `xq` is itself sorted. Queries at or below `x[1]` clamp to `y[1]` and queries at or above `x[end]` clamp to `y[end]`. A degenerate bracket where `x1 == x0` sets the interpolation weight to `0.0` instead of dividing by zero.

## Theory & Math

Inside a bracket $[x_j, x_{j+1}]$ the returned value is the convex combination

$$y(q) = y_j + w\,(y_{j+1} - y_j), \qquad w = \frac{q - x_j}{x_{j+1} - x_j}$$

where $q$ is the query abscissa, $x_j$ and $x_{j+1}$ are the bracketing knots and $y_j$, $y_{j+1}$ the corresponding ordinates. Outside the knot range $y(q)$ is held constant at the nearest endpoint value.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Vector{Float64} | n/a | yes | Positional argument `x`. |
| in | `y` | Vector{Float64} | n/a | yes | Positional argument `y`. |
| in | `xq` | Vector{Float64} | n/a | yes | Positional argument `xq`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_interp_linear`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.comparison_metrics__compare_time_series|_compare_time_series]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:268-268`
- [[analysis.comparison_metrics__rates|_rates]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:77-77`
- [[envana.ana_comparison_metrics_compare_orbit_curve|_compare_orbit_curve]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:186-186`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The forward-only cursor `j` is never reset between query points, so an unsorted `xq` will be interpolated against stale brackets and produce wrong values without any error. The abscissa `x` is assumed strictly increasing; duplicate or decreasing knots are not detected beyond the equal-bracket guard. Extrapolation is clamped rather than linear, which masks the fact that the simulation did not cover the requested span. Arguments are restricted to concrete `Vector{Float64}`, and `@inbounds` disables bounds checking in the hot loop.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/comparison_metrics.jl` line 1.
