---
id: analysis.manifest_parsing__parse_truth_mask
label: _parse_truth_mask
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_truth_mask
  lines:
  - 197
  - 197
inputs:
- id: raw
  type: String
  units: n/a
  required: true
  description: Positional argument `raw`.
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
  type: Symbol
  units: n/a
  description: Return value of `_parse_truth_mask`.
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

# _parse_truth_mask

## Purpose
Maps the truth-mask string onto `:none`, `:nightside` or `:dayside`, selecting which portion of each orbit the telemetry comparison scores.

## Design & Implementation
Lowercases and strips the input, treats an empty string as `:none`, accepts the three names exactly, and raises `ArgumentError` listing the accepted values otherwise. `@inline` with a `::Symbol` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_parse_truth_mask`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:701-701`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Unlike the other vocabulary parsers in this file it accepts no aliases, so `night` or `day` is rejected; the mask semantics themselves are implemented downstream in the telemetry loader.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 197.
