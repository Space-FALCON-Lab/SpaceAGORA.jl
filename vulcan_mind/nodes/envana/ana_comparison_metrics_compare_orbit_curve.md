---
id: envana.ana_comparison_metrics_compare_orbit_curve
label: _compare_orbit_curve
kind: function
source:
  file: src/analysis/verification/telemetry_verification/comparison_metrics.jl
  symbol: _compare_orbit_curve
  lines:
  - 97
  - 225
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying interpolation, decay diagnostics,
    and statistics.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: summary_row
  type: NamedTuple
  units: km
  description: Accuracy summary carrying mae, rmse, max_abs, p95_abs, bias, nmae,
    nrmse and coverage for one event.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- envana
origin: agent
---
# _compare_orbit_curve

## Purpose
`_compare_orbit_curve` scores one simulated per-orbit curve against the matching telemetry curve for a named scenario and event, returning both a summary named tuple and the per-point error rows used to build the error table.

## Theory & Math
With `n` compared points, simulated values `s_j` in kilometres, telemetry values `t_j` in kilometres, and residuals `e_j = s_j + b - t_j` where `b` is the calibration bias, the routine reports `mae = mean(|e|)`, `rmse = sqrt(mean(e^2))`, `max_abs = max(|e|)`, the 95th percentile of `|e|`, and the signed mean bias, all in kilometres. It also reports the scale-free forms `nmae = mean(|e|) / R` and `nrmse = rmse / R`, where `R` is the peak-to-peak range of the telemetry series in kilometres, making the two normalised metrics dimensionless and comparable across events whose altitudes differ by orders of magnitude.

## Model & Assumptions
Telemetry and simulation may sit on different orbit-number axes, so the simulated series is resampled onto the telemetry axis by `_interp_linear` before differencing. When `mask_to_sim_span` is set and a `sim_axis` is supplied, only telemetry points inside `[sim_axis[1] - 0.5, sim_axis[end] + 0.5]` are scored, the half-orbit margin giving grace at each end; `coverage` then reports the masked fraction of the full telemetry series. Without masking, telemetry beyond the simulated span is compared against the clamped final simulated value, which the source describes as accumulation bookkeeping rather than trajectory error.

## Design & Implementation
Two guard branches short-circuit degenerate input: an empty telemetry series raises immediately, and a zero-length simulated series returns a named tuple with every metric set to `Inf` and `coverage` zero, merged with `_DECAY_DIAGNOSTIC_EMPTY` so the field set stays uniform. Apoapsis curves with at least three samples additionally receive `_apo_decay_diagnostic` fields. Returning a `NamedTuple` keeps the field names available to `DataFrame` construction without a separate schema.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying interpolation, decay diagnostics, and statistics. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `summary_row` | NamedTuple | km | — | Accuracy summary carrying mae, rmse, max_abs, p95_abs, bias, nmae, nrmse and coverage for one event. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:97-97`
- [[envana.ana_error_tables_orbit_rows_errors|_orbit_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:27-27`

**Downstream**

- `callees` → [[analysis.comparison_metrics__interp_linear|_interp_linear]] · `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:186-186`
- `callees` → [[analysis.comparison_metrics__normalized_axis|_normalized_axis]] · `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:184-184`
<!-- vulcan:connections:end -->

## Limitations
Linear interpolation between orbit samples underestimates curvature when the telemetry cadence is coarse. Normalising by the telemetry peak-to-peak range makes `nmae` and `nrmse` unstable for nearly flat series, where `R` approaches zero. The `Inf` sentinel for empty simulations propagates into any aggregate statistic computed over multiple events.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/comparison_metrics.jl:97-225`, including the masking branch, the degenerate-input returns, and the metric expressions at lines 206 through 211.
