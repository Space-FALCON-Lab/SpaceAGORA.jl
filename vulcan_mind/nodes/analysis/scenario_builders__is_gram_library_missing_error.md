---
id: analysis.scenario_builders__is_gram_library_missing_error
label: _is_gram_library_missing_error
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _is_gram_library_missing_error
  lines:
  - 289
  - 289
inputs:
- id: err
  type: Any
  units: n/a
  required: true
  description: Positional argument `err`.
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
  description: Return value of `_is_gram_library_missing_error`.
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

# _is_gram_library_missing_error

## Purpose
Recognises the specific failure of GRAM being unavailable, so the builder can fall back instead of rethrowing.

## Design & Implementation
Renders the exception with `showerror`, lowercases it, and checks for the substrings `gram shared library` or `gram_lib`. `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `err` | Any | n/a | yes | Positional argument `err`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_is_gram_library_missing_error`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_required_gram_density_model|_make_required_gram_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:339-339`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Matching on message text is fragile; a rewording of the GRAMSuite error would defeat the fallback and turn a recoverable condition into a hard failure.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 289.
