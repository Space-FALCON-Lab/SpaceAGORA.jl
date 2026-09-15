---
id: analysis.reporting__display_value_scale
label: _display_value_scale
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _display_value_scale
  lines:
  - 40
  - 40
inputs:
- id: value_units
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `value_units`.
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
  description: Return value of `_display_value_scale`.
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

# _display_value_scale

## Purpose
Returns the multiplicative factor that converts an internal value to its display unit, pairing with `_display_value_units`: 1000.0 for `km/s` to `m/s`, 1.0 for everything else.

## Design & Implementation
The input `value_units::AbstractString` is normalised with `lowercase(strip(String(...)))` and compared to `"km/s"`. Returns `Float64`: `1e3` on match, `1.0` otherwise. Called once per row to build a `scales` vector that is then broadcast across every numeric column selected for display.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value_units` | AbstractString | n/a | yes | Positional argument `value_units`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_display_value_scale`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.reporting__append_display_error_columns_bang|_append_display_error_columns!]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:85-85`
- [[analysis.reporting__append_display_metric_columns_bang|_append_display_metric_columns!]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:50-50`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A mismatch between this function and `_display_value_units` would silently produce wrong numbers with the wrong label; nothing enforces consistency other than both using the same literal. Because `km` maps to scale 1.0, `*_display` columns for distance events are numerically identical to the `*_km` columns.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 40.
