---
id: analysis.manifest_parsing__require_key
label: _require_key
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _require_key
  lines:
  - 56
  - 56
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
  description: Return value of `_require_key`. Returns `tbl[key]`.
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

# _require_key

## Purpose
The base accessor for mandatory manifest fields, producing an error message that names both the key and its location.

## Design & Implementation
Returns `tbl[key]` if present, otherwise throws `ArgumentError("Missing 'key' in context")`. Every `_require_*` helper builds on it. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_require_key`. Returns `tbl[key]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:524-524`
- [[analysis.manifest_parsing__parse_vec3|_parse_vec3]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:452-452`
- [[analysis.manifest_parsing__require_float|_require_float]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:72-72`
- [[analysis.manifest_parsing__require_int|_require_int]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:76-76`
- [[analysis.manifest_parsing__require_str|_require_str]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:68-68`
- [[analysis.manifest_parsing__require_str_vector|_require_str_vector]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:110-110`
- [[analysis.manifest_parsing__require_table|_require_table]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:62-62`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It returns the raw TOML value, so callers must still convert and validate type.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 56.
