---
id: analysis.reporting__display_value_units
label: _display_value_units
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _display_value_units
  lines:
  - 35
  - 35
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
  type: String
  units: n/a
  description: Return value of `_display_value_units`.
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

# _display_value_units

## Purpose
Maps an internal unit string to the unit shown in human-facing report columns, converting `km/s` to `m/s` so small velocity errors are readable, and leaving every other unit untouched.

## Design & Implementation
Takes any `AbstractString`, normalises it via `lowercase(strip(String(value_units)))` into `token`, and returns `"m/s"` when `token == "km/s"`, otherwise the original (un-normalised) string. The companion `_display_value_scale` must agree with this mapping because the two are always applied together by `_append_display_metric_columns!` and `_append_display_error_columns!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `value_units` | AbstractString | n/a | yes | Positional argument `value_units`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_display_value_units`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.reporting__append_display_error_columns_bang|_append_display_error_columns!]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:85-85`
- [[analysis.reporting__append_display_metric_columns_bang|_append_display_metric_columns!]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:51-51`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only one conversion is recognised; `km` is not converted to `m` even though the metric columns are named `*_km`, so distance errors stay in kilometres while velocity errors switch units. Variants like `"km / s"` or `"kps"` are not matched. The mapping and its scale are hard-coded in two separate functions with no shared table.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 35.
