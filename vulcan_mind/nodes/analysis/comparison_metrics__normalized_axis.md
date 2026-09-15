---
id: analysis.comparison_metrics__normalized_axis
label: _normalized_axis
kind: function
source:
  file: src/analysis/verification/telemetry_verification/comparison_metrics.jl
  symbol: _normalized_axis
  lines:
  - 24
  - 24
inputs:
- id: n
  type: Int
  units: n/a
  required: true
  description: Positional argument `n`.
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
  description: Return value of `_normalized_axis`. Returns `collect(range(0.0, 1.0,
    length=n))`.
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

# _normalized_axis

## Purpose

Produces a unitless progress axis spanning zero to one with `n` evenly spaced samples. It is used when a telemetry series and a simulated series have no common physical abscissa, so both are placed on a shared normalised index before interpolation and differencing.

## Design & Implementation

`_normalized_axis(n::Int)` returns the single-element vector `[0.0]` when `n <= 1`, avoiding the zero-length step that `range` would otherwise produce, and otherwise returns `collect(range(0.0, 1.0, length=n))`. It is marked `@inline`. `_compare_orbit_curve` calls it twice, once for the telemetry length and once for the simulation length, whenever no explicit `sim_axis` is supplied.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n` | Int | n/a | yes | Positional argument `n`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_normalized_axis`. Returns `collect(range(0.0, 1.0, length=n))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_comparison_metrics_compare_orbit_curve|_compare_orbit_curve]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:184-184`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Normalising both series to the unit interval implicitly assumes they cover the same physical span, so a simulation that stops early is stretched to match the full telemetry duration and the resulting errors understate the true divergence. The masking and axis-aware branches exist precisely to avoid that, making this a fallback rather than the preferred comparison. A negative or zero `n` silently yields `[0.0]` instead of an empty result, and the returned vector is heap-allocated on every call.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/comparison_metrics.jl` line 24.
