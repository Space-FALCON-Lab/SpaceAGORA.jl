---
id: analysis.manifest_parsing__optional_float_tuple
label: _optional_float_tuple
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _optional_float_tuple
  lines:
  - 286
  - 286
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
- id: n
  type: Int
  units: n/a
  required: true
  description: Positional argument `n`.
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
  type: Nothing
  units: n/a
  description: Return value of `_optional_float_tuple`. Returns `nothing` or `values`.
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

# _optional_float_tuple

## Purpose
Fetches an optional fixed-length numeric array as an `NTuple`, used for state overrides, GRAM scales and Mars parameters.

## Design & Implementation
Returns `nothing` if absent. Otherwise requires an array of exactly `n` elements, raising with the expected and actual lengths, and builds the tuple with `ntuple` and `Float64` conversion. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `n` | Int | n/a | yes | Positional argument `n`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_optional_float_tuple`. Returns `nothing` or `values`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:602-602`
- [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:340-340`
- [[analysis.manifest_parsing__parse_ic_offset|_parse_ic_offset]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:191-191`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:293-293`
<!-- vulcan:connections:end -->

## Limitations
The return type depends on the runtime `n`, so this is not type-stable; acceptable for parse-time code but not for a hot path.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 286.
