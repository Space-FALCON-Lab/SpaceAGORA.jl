---
id: analysis.manifest_parsing__parse_positive_int_env
label: _parse_positive_int_env
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_positive_int_env
  lines:
  - 11
  - 11
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Int
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
  type: Int
  units: n/a
  description: Return value of `_parse_positive_int_env`.
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

# _parse_positive_int_env

## Purpose
Reads a positive integer from an environment variable, falling back to a default when unset and failing loudly when malformed.

## Design & Implementation
Strips the variable, returns `default` when empty, parses with `parse(Int, ...)` inside a `try` that converts a parse failure into an `ArgumentError` naming the variable and raw value, and rejects non-positive results with a second `ArgumentError`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Int | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_parse_positive_int_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__telemetry_solver_maxiters|_telemetry_solver_maxiters]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:37-37`
- [[analysis.manifest_parsing__telemetry_solver_retry_maxiters|_telemetry_solver_retry_maxiters]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:42-42`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Unlike `_safe_parse_bool`, this one throws on bad input, which is the right choice for iteration caps but means the two environment helpers have opposite error philosophies.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 11.
