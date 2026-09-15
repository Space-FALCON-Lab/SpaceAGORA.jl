---
id: analysis.manifest_parsing__require_str_vector
label: _require_str_vector
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _require_str_vector
  lines:
  - 109
  - 109
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
  type: Vector{String}
  units: n/a
  description: Return value of `_require_str_vector`.
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

# _require_str_vector

## Purpose
Fetches a mandatory array of strings, used for the scenario's list of compared event names.

## Design & Implementation
Calls `_require_key`, requires the value to be an `AbstractVector` with an error naming the key and context, and maps `String` over the elements into a `Vector{String}`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{String} | n/a | — | Return value of `_require_str_vector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:536-536`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__require_key|_require_key]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:110-110`
<!-- vulcan:connections:end -->

## Limitations
An empty array passes validation, which for the events list means no channels are compared and the scenario trivially succeeds without any warning that nothing was checked.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 109.
