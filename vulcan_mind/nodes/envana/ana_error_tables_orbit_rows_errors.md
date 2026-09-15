---
id: envana.ana_error_tables_orbit_rows_errors
label: _orbit_rows_errors
kind: function
source:
  file: src/analysis/verification/telemetry_verification/error_tables.jl
  symbol: _orbit_rows_errors
  lines:
  - 1
  - 59
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying curve comparison and apoapsis
    decay diagnostics.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: rows_errors
  type: Tuple{Vector,Vector}
  units: km
  description: Periapsis and apoapsis summary rows paired with their per-point error
    tables.
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
# _orbit_rows_errors

## Purpose
`_orbit_rows_errors` turns one orbit-events scenario into the two comparison records the study reports: a periapsis altitude record and an apoapsis altitude record, each with a summary row and a per-point error table.

## Theory & Math
Both records compare simulated extrema altitudes against telemetry extrema on an orbit-number axis. When the scenario does not carry an explicit epoch offset, the simulated axis is reconstructed as `x_j = x_1 + j * delta`, where `x_1` is the first telemetry orbit number and `delta = median(diff(orbit))` is the robust orbit-number spacing; the median is used so a single dropped orbit does not stretch the axis. Otherwise the axis starts at `cfg.epoch_orbit_offset` with unit spacing. The apoapsis record additionally receives the secular decay diagnostic, whose slope is expressed in metres per day.

## Model & Assumptions
The routine assumes telemetry supplies both periapsis and apoapsis series and that the simulated extrema arrays are indexable in orbit order. Per-event calibration biases `peri_bias` and `apo_bias` are passed through to the comparison so that a constant altitude offset is removed before scoring. The `mask_to_sim` flag decides whether telemetry beyond the simulated span is excluded. The decay diagnostic is attached only when at least three apoapsis samples exist, since the harmonic fit needs an overdetermined system.

## Design & Implementation
Two calls to `_compare_orbit_curve` produce the periapsis and apoapsis results with event labels `"peri"` and `"apo"`. Burn orbit numbers come from `cfg.maneuver_orbit_numbers_campaign` when it is non-empty and from `cfg.maneuver_orbit_numbers` otherwise, so a campaign definition takes precedence over the plain list. The apoapsis summary is extended with `merge(apo_summary, _apo_decay_diagnostic(...))`, which keeps the named tuple immutable while adding fields. The bias is added to the simulated apoapsis altitudes before the diagnostic runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying curve comparison and apoapsis decay diagnostics. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `rows_errors` | Tuple{Vector,Vector} | km | — | Periapsis and apoapsis summary rows paired with their per-point error tables. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:166-166`

**Downstream**

- `callees` → [[analysis.comparison_metrics__apo_decay_diagnostic|_apo_decay_diagnostic]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:50-50`
- `callees` → [[analysis.telemetry_loading__extract_extrema_series|_extract_extrema_series]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:8-8`
- `callees` → [[analysis.telemetry_loading__load_telemetry_curve|_load_telemetry_curve]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:9-9`
- `callees` → [[envana.ana_comparison_metrics_compare_orbit_curve|_compare_orbit_curve]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:27-27`
<!-- vulcan:connections:end -->

## Limitations
Only periapsis and apoapsis are tabulated; intermediate along-track behaviour is invisible to this path. Reconstructing the simulated axis from the telemetry cadence assumes the simulation produced one extremum per telemetry orbit, which fails if the simulation terminates early or skips an apsis. Fewer than two telemetry points force the fallback spacing of 1.0, which silently mislabels the axis.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/error_tables.jl:1-59`, including the axis reconstruction branch and the apoapsis decay merge.
