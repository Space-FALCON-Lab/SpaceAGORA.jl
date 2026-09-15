---
id: analysis.ic_fit__ic_fit_series
label: _ic_fit_series
kind: function
source:
  file: src/analysis/verification/telemetry_verification/ic_fit.jl
  symbol: _ic_fit_series
  lines:
  - 16
  - 16
inputs:
- id: errors_csv
  type: String
  units: n/a
  required: true
  description: Positional argument `errors_csv`.
- id: scenario_name
  type: String
  units: n/a
  required: true
  description: Positional argument `scenario_name`.
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
  description: Return value of `_ic_fit_series`. Returns `out`.
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

# _ic_fit_series

## Purpose

Reads a verification errors CSV produced by a fit run and extracts, per position event, the per-sample simulated value and residual keyed by telemetry sample time, giving the differential-correction solver the raw numbers it differences.

## Design & Implementation

Loads `errors_csv` with `CSV.read(..., DataFrame)`, then filters to rows whose `scenario` column equals `scenario_name` and whose `event` is one of `_IC_FIT_EVENTS` = `("state_x_time", "state_y_time", "state_z_time")`, throwing `ArgumentError` when nothing survives. For each of the three events it builds a `Dict{Float64, NTuple{2, Float64}}` mapping `telemetry_axis` (the comparison sample abscissa) to `(sim_interp_value_km, error_km)`, and returns those three dictionaries in a `Dict{String, ...}` keyed by event name. Dictionary keying, rather than positional vectors, is what lets the caller intersect sample sets across the baseline and the six perturbed runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `errors_csv` | String | n/a | yes | Positional argument `errors_csv`. |
| in | `scenario_name` | String | n/a | yes | Positional argument `scenario_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_ic_fit_series`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.ic_fit__ic_fit_run|_ic_fit_run]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:61-61`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Sample abscissae are used as raw `Float64` dictionary keys, so two runs whose comparison times differ by a floating-point ulp will not be recognised as the same sample and will silently drop out of the common set. Duplicate rows for one event and abscissa are collapsed by the comprehension with last-write-wins and no warning. The whole CSV is read into memory before filtering, and the three required column names (`telemetry_axis`, `sim_interp_value_km`, `error_km`) are assumed present.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/ic_fit.jl` line 16.
