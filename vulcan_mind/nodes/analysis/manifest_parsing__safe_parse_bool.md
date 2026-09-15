---
id: analysis.manifest_parsing__safe_parse_bool
label: _safe_parse_bool
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _safe_parse_bool
  lines:
  - 1
  - 1
inputs:
- id: raw
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `raw`.
- id: default
  type: Bool
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
  type: Bool
  units: n/a
  description: Return value of `_safe_parse_bool`.
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

# _safe_parse_bool

## Purpose
Interprets a boolean-like string from the environment or command line, returning a default rather than failing on unrecognised text.

## Design & Implementation
Lowercases and strips the input, returns true for `1`, `true`, `yes` or `on`, false for `0`, `false`, `no` or `off`, and `default` for anything else. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | AbstractString | n/a | yes | Positional argument `raw`. |
| in | `default` | Bool | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_safe_parse_bool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__libraryless_gram_surrogate_enabled|_libraryless_gram_surrogate_enabled]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:295-295`
- [[analysis.types_verificationrequest|VerificationRequest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:262-262`
- [[envana.ana_manifest_parsing_parse_cli|parse_cli]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:717-717`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A misspelling silently yields the default, so a typo in a CI variable can flip an enforcement flag without any diagnostic.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 1.
