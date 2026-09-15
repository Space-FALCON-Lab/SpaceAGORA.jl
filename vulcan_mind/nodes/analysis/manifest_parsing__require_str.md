---
id: analysis.manifest_parsing__require_str
label: _require_str
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _require_str
  lines:
  - 67
  - 67
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
  type: String
  units: n/a
  description: Return value of `_require_str`.
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

# _require_str

## Purpose
Fetches a mandatory string field from a manifest table, used for scenario names, planet names, file paths and enumerated mode strings before they are parsed further.

## Design & Implementation
Wraps `_require_key` in `String(...)` and is declared `@inline` with a `::String` return. The missing-key case therefore produces the same `Missing 'key' in context` message as every other required accessor, keeping manifest errors uniform.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_require_str`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:533-533`
- [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:303-303`
- [[analysis.manifest_parsing__parse_units|_parse_units]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:490-490`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__require_key|_require_key]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:68-68`
<!-- vulcan:connections:end -->

## Limitations
`String` of a non-string TOML value such as a number throws a `MethodError` rather than an `ArgumentError`, and that failure does not name the key, unlike the missing-key message.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 67.
