---
id: analysis.manifest_parsing__optional_float
label: _optional_float
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _optional_float
  lines:
  - 83
  - 83
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
- id: default
  type: Float64
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: Return value of `_optional_float`.
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

# _optional_float

## Purpose
Fetches a float field with a default when absent, used widely for tolerances, scale factors, thrust, Isp and effector rates.

## Design & Implementation
Returns `Float64(tbl[key])` if present and `default` otherwise, declared `@inline` with a `::Float64` return. Because conversion goes through `Float64`, integer literals in the manifest are accepted for float-valued keys.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `default` | Float64 | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_optional_float`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:551-551`
- [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:308-308`
- [[analysis.manifest_parsing__parse_calibration_config|_parse_calibration_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:402-402`
- [[analysis.manifest_parsing__parse_event_tolerance|_parse_event_tolerance]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:503-503`
- [[analysis.manifest_parsing__parse_maneuver_config|_parse_maneuver_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:279-279`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations
No validation of sign or finiteness; callers such as the manoeuvre and calibration parsers add their own positivity and ordering checks after reading.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 83.
