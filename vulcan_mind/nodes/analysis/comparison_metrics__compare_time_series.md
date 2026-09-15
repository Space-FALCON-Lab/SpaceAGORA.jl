---
id: analysis.comparison_metrics__compare_time_series
label: _compare_time_series
kind: function
source:
  file: src/analysis/verification/telemetry_verification/comparison_metrics.jl
  symbol: _compare_time_series
  lines:
  - 227
  - 227
inputs:
- id: scenario
  type: String
  units: n/a
  required: true
  description: Positional argument `scenario`.
- id: event
  type: String
  units: n/a
  required: true
  description: Positional argument `event`.
- id: telemetry_time
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `telemetry_time`.
- id: telemetry_values
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `telemetry_values`.
- id: sim_time
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `sim_time`.
- id: sim_values
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `sim_values`.
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
  description: Return value of `_compare_time_series`. Returns `(`.
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

# _compare_time_series

## Purpose

Scores a simulated time-history against a flight telemetry time-history for one scenario and event, returning both a summary metric row and a per-sample residual table that downstream reporting writes out for plotting.

## Design & Implementation

`_compare_time_series(scenario, event, telemetry_time, telemetry_values, sim_time, sim_values)` errors immediately when the telemetry series is empty. If the simulation produced no samples it returns an all-`Inf` metric row with `coverage = 0.0`, the `_DECAY_DIAGNOSTIC_EMPTY` placeholder fields and an empty `DataFrame` with the standard column schema, so every row shares one shape. Otherwise the simulation is resampled onto `telemetry_time` with `_interp_linear`, residuals are `sim_interp .- telemetry_values`, and the summary reports `mae_km`, `rmse_km`, `max_abs_km`, `p95_abs_km` (the 0.95 quantile of absolute error) and `bias_km` (the signed mean). Normalised forms `nmae` and `nrmse` divide by the telemetry peak-to-peak range floored at `1e-9`. The returned `DataFrame` carries the axis, both values, the interpolated simulation twice (raw and biased are identical here) and the signed error.

## Theory & Math

With residuals $e_i = \hat{y}(t_i) - y_i$ over the $N$ telemetry samples, where $\hat{y}$ is the interpolated simulation and $y_i$ the telemetry value in km,

$$\text{MAE} = \frac{1}{N}\sum_{i=1}^{N} |e_i|, \qquad \text{RMSE} = \sqrt{\frac{1}{N}\sum_{i=1}^{N} e_i^2}, \qquad b = \frac{1}{N}\sum_{i=1}^{N} e_i$$

and the normalised variants divide by $R = \max(\max_i y_i - \min_i y_i,\; 10^{-9})$, making them dimensionless.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `scenario` | String | n/a | yes | Positional argument `scenario`. |
| in | `event` | String | n/a | yes | Positional argument `event`. |
| in | `telemetry_time` | Vector{Float64} | n/a | yes | Positional argument `telemetry_time`. |
| in | `telemetry_values` | Vector{Float64} | n/a | yes | Positional argument `telemetry_values`. |
| in | `sim_time` | Vector{Float64} | n/a | yes | Positional argument `sim_time`. |
| in | `sim_values` | Vector{Float64} | n/a | yes | Positional argument `sim_values`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_compare_time_series`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:179-179`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl`

**Downstream**

- `callees` → [[analysis.comparison_metrics__interp_linear|_interp_linear]] · `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:268-268`
<!-- vulcan:connections:end -->

## Limitations

`coverage` is computed as `min(n_tel, n_sim) / n_tel`, a sample-count ratio rather than a true temporal-overlap fraction, so a short simulation that happens to be densely sampled reports full coverage while its extrapolated tail is clamped to the final simulated value by `_interp_linear`. That clamping means residuals beyond the simulated span measure bookkeeping rather than trajectory error, and unlike `_compare_orbit_curve` this function offers no span masking or bias option. Column names are fixed to kilometre units regardless of the physical quantity being compared, and the empty-simulation branch reads `telemetry_time[1]` and `telemetry_time[end]`, which is safe only because the empty-telemetry case errors first.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/comparison_metrics.jl` line 227.
