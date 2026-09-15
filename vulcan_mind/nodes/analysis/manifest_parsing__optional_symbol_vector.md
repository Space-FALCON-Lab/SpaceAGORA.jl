---
id: analysis.manifest_parsing__optional_symbol_vector
label: _optional_symbol_vector
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _optional_symbol_vector
  lines:
  - 95
  - 95
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
  type: Vector{Symbol}
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
  type: Vector{Symbol}
  units: n/a
  description: Return value of `_optional_symbol_vector`.
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

# _optional_symbol_vector

## Purpose
Fetches an optional array of names as lowercase symbols, used for the calibration profile list.

## Design & Implementation
Returns a copy of `default` if absent. Otherwise requires an array, converts each element to a stripped lowercase `Symbol`, and rejects an empty result with an `ArgumentError` that specifically says calibration profiles cannot be empty. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `default` | Vector{Symbol} | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Symbol} | n/a | — | Return value of `_optional_symbol_vector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_calibration_config|_parse_calibration_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:393-393`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:103-103`
<!-- vulcan:connections:end -->

## Limitations
The error message is calibration-specific even though the function is generic, so reusing it for another key would produce a misleading message.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 95.
