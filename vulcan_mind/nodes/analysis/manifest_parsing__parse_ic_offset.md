---
id: analysis.manifest_parsing__parse_ic_offset
label: _parse_ic_offset
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_ic_offset
  lines:
  - 190
  - 190
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
  type: NTuple{3,
  units: n/a
  description: Return value of `_parse_ic_offset`.
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

# _parse_ic_offset

## Purpose
Parses an optional three-component initial-condition offset, defaulting to zero and rejecting non-finite entries.

## Design & Implementation
Uses `_optional_float_tuple` with `n = 3`, returns `(0.0, 0.0, 0.0)` when absent, and raises if any component is non-finite. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NTuple{3, | n/a | — | Return value of `_parse_ic_offset`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:699-699`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__optional_float_tuple|_optional_float_tuple]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:191-191`
<!-- vulcan:connections:end -->

## Limitations
Units are implied by the key name (`_m` or `_mps`) rather than checked.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 190.
