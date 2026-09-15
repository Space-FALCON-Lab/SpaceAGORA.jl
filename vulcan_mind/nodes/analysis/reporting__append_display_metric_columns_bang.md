---
id: analysis.reporting__append_display_metric_columns_bang
label: _append_display_metric_columns!
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _append_display_metric_columns!
  lines:
  - 45
  - 45
inputs:
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
  description: Return value of `_append_display_metric_columns!`; mutates `summary_df`
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

# _append_display_metric_columns!

## Purpose
Augments the per-event summary DataFrame with a `value_units_display` column and `*_display` versions of every error-metric and tolerance-limit column, so downstream CSV and plots can present velocity errors in m/s while keeping the canonical `*_km` columns intact.

## Design & Implementation
Mutates `summary_df::DataFrame` in place and returns `nothing`. Early-returns when `nrow == 0` or when the `value_units` column is absent. It materialises `units_raw` as `String`s, computes `scales` via `_display_value_scale` and `units_display` via `_display_value_units`, assigns `summary_df.value_units_display`, then invokes the nested `_add_scaled_column!` closure seven times for `mae_km`, `rmse_km`, `max_abs_km`, `p95_abs_km`, `bias_km`, `limit_max_abs_km`, `limit_max_rmse_km` mapping each to its `_display` counterpart. Columns that do not exist are skipped silently.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `summary_df` | DataFrame | n/a | yes | Positional argument `summary_df`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_append_display_metric_columns!`; mutates `summary_df` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:431-431`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- `callees` → [[analysis.reporting__display_value_scale|_display_value_scale]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:50-50`
- `callees` → [[analysis.reporting__display_value_units|_display_value_units]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
`nmae` and `limit_nmae` are dimensionless and deliberately not scaled, but nothing documents that in the frame. Any `missing` in a scaled column raises a conversion error. Existing `*_display` columns are overwritten without warning. The function assumes `value_units` entries are convertible with `String(v)`.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 45.
