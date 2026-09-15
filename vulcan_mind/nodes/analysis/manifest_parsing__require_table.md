---
id: analysis.manifest_parsing__require_table
label: _require_table
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _require_table
  lines:
  - 61
  - 61
inputs:
- id: tbl
  type: Any
  units: n/a
  required: true
  description: Positional argument `tbl`.
- id: key
  type: String
  units: n/a
  required: true
  description: Positional argument `key`.
- id: context
  type: String
  units: n/a
  required: true
  description: Positional argument `context`.
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
  description: Return value of `_require_table`. Returns `value`.
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

# _require_table

## Purpose
Fetches a mandatory sub-table from a manifest table, rejecting a scalar or array at that key.

## Design & Implementation
Calls `_require_key` then checks `isa AbstractDict`, throwing `ArgumentError("Expected table ...")` otherwise. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_require_table`. Returns `value`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:541-541`
- [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:301-301`
- [[analysis.manifest_parsing__parse_calibration_config|_parse_calibration_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:392-392`
- [[analysis.manifest_parsing__parse_event_tolerance|_parse_event_tolerance]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:499-499`
- [[analysis.manifest_parsing__parse_maneuver_config|_parse_maneuver_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:219-219`
- [[analysis.manifest_parsing__parse_spacecraft_config|_parse_spacecraft_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:472-472`
- [[analysis.manifest_parsing__parse_tolerances|_parse_tolerances]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:508-508`
- [[analysis.manifest_parsing__parse_units|_parse_units]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:489-489`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__require_key|_require_key]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:62-62`
<!-- vulcan:connections:end -->

## Limitations
An empty table passes, so a section that is present but blank fails later at its first required field rather than here.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 61.
