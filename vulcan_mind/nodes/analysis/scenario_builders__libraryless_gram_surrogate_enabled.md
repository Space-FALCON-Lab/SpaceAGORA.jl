---
id: analysis.scenario_builders__libraryless_gram_surrogate_enabled
label: _libraryless_gram_surrogate_enabled
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _libraryless_gram_surrogate_enabled
  lines:
  - 294
  - 294
inputs:
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
  description: Return value of `_libraryless_gram_surrogate_enabled`.
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

# _libraryless_gram_surrogate_enabled

## Purpose
Reads the switch that allows a telemetry run to substitute the offline GRAM surrogate when the native library is missing.

## Design & Implementation
Parses `SPACEAGORA_TELEMETRY_ALLOW_GRAM_OFFLINE_NO_LIB` through `_safe_parse_bool` with a default of true, so the fallback is on unless explicitly disabled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_libraryless_gram_surrogate_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__try_libraryless_gram_surrogate|_try_libraryless_gram_surrogate]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:303-303`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__safe_parse_bool|_safe_parse_bool]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:295-295`
<!-- vulcan:connections:end -->

## Limitations
Defaulting to enabled means a CI machine without the library quietly produces surrogate-based results; a benchmark that must use native GRAM has to set this to false or set `gram_offline_surrogate` to `off`.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 294.
