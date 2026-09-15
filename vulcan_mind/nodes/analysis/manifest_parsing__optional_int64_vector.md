---
id: analysis.manifest_parsing__optional_int64_vector
label: _optional_int64_vector
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _optional_int64_vector
  lines:
  - 124
  - 124
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
  type: Vector{Int64}
  units: n/a
  description: Return value of `_optional_int64_vector`.
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

# _optional_int64_vector

## Purpose
Fetches an optional array of integers as `Vector{Int64}`, used for the manoeuvre orbit-number list.

## Design & Implementation
Returns an empty `Int64` vector if the key is absent; otherwise requires an `AbstractVector` and converts each element through `Int64`, so TOML integers of any width are accepted. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Int64} | n/a | — | Return value of `_optional_int64_vector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_maneuver_config|_parse_maneuver_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:220-220`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Carries the same fixed `manifest scenario table` context string in its error as the string-vector sibling, and does not check that entries are positive — the manoeuvre parser does that afterwards.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 124.
