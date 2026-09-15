---
id: analysis.reporting__append_display_error_columns_bang
label: _append_display_error_columns!
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _append_display_error_columns!
  lines:
  - 75
  - 75
inputs:
- id: errors_df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `errors_df`.
- id: summary_df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `summary_df`.
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
  type: Nothing
  units: n/a
  description: Return value of `_append_display_error_columns!`; mutates `errors_df`
    in place.
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

# _append_display_error_columns!

## Purpose
Adds display-unit columns to the per-sample error DataFrame by looking up each row's `(scenario, event)` pair in the summary frame to find the unit and scale that apply, so the sample-level errors match the aggregate metrics' presentation.

## Design & Implementation
Mutates `errors_df` and returns `nothing`; `summary_df` is read-only. Guards return early when `errors_df` is empty or lacks `scenario`/`event` columns. A `Dict{Tuple{String,String},Tuple{Float64,String}}` named `unit_map` is populated from `summary_df` rows (only if it has `scenario`, `event`, `value_units`), storing `(_display_value_scale(u), _display_value_units(u))`. A second `@inbounds` pass over `1:nrow(errors_df)` fills `scales` and `units_display`, defaulting to `(1.0, "km")` when the key is missing. After assigning `errors_df.value_units_display`, the nested `_add_scaled_column!` is applied to `telemetry_value_km`, `sim_interp_value_km`, and `error_km`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `errors_df` | DataFrame | n/a | yes | Positional argument `errors_df`. |
| in | `summary_df` | DataFrame | n/a | yes | Positional argument `summary_df`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_append_display_error_columns!`; mutates `errors_df` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:432-432`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:105-105`
- `callees` → [[analysis.reporting__add_scaled_column_bang|_add_scaled_column!]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:100-100`
- `callees` → [[analysis.reporting__display_value_scale|_display_value_scale]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:85-85`
- `callees` → [[analysis.reporting__display_value_units|_display_value_units]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:85-85`
<!-- vulcan:connections:end -->

## Limitations
Rows whose `(scenario, event)` pair is absent from the summary silently get km with scale 1.0, which mislabels velocity events if the summary was filtered. Duplicate summary rows for the same key let the last one win. The `String(...)` conversions throw on `missing` scenario or event cells.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 75.
