---
id: analysis.telemetry_loading__differentiate_series
label: _differentiate_series
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _differentiate_series
  lines:
  - 411
  - 411
inputs:
- id: values
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `values`.
- id: time_s
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `time_s`.
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
  type: Vector{Float64}
  units: n/a
  description: Return value of `_differentiate_series`.
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

# _differentiate_series

## Purpose
Numerically differentiates a sampled scalar series with respect to non-uniform time stamps, producing velocity components (km/s) from telemetry position columns when no measured velocity is available.

## Theory & Math
$$\dot v_1 = \frac{v_2 - v_1}{t_2 - t_1},\qquad \dot v_i = \frac{v_{i+1} - v_{i-1}}{t_{i+1} - t_{i-1}}\ (1<i<n),\qquad \dot v_n = \frac{v_n - v_{n-1}}{t_n - t_{n-1}}$$ with each denominator replaced by $\max(\Delta t, \epsilon_{mach})$.

## Design & Implementation
Checks `length(values) == length(time_s)` and `n >= 2`, throwing `ArgumentError` otherwise. Allocates `dv::Vector{Float64}(undef, n)` and, inside `@inbounds`, uses a forward difference at index 1, a central difference `(v[i+1] - v[i-1]) / (t[i+1] - t[i-1])` for interior points, and a backward difference at index `n`. Every time step is floored at `eps(Float64)` with `max(dt, eps)` so repeated time stamps produce very large finite values rather than `Inf` or `NaN`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `values` | Vector{Float64} | n/a | yes | Positional argument `values`. |
| in | `time_s` | Vector{Float64} | n/a | yes | Positional argument `time_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Float64} | n/a | — | Return value of `_differentiate_series`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.telemetry_loading__load_time_aligned_telemetry|_load_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:373-373`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The interior formula is only first-order accurate on non-uniform grids and second-order only when spacing is uniform. Flooring `dt` at `eps` turns duplicate time stamps into derivatives of order `1e16`, which then contaminate any downstream statistic instead of raising an error. No smoothing is applied, so measurement noise in position is amplified by `1/dt`. Time stamps are assumed sorted ascending; a descending pair yields `eps` in the denominator, not a negative step.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 411.
