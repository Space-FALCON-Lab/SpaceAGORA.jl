---
id: analysis.manifest_parsing__telemetry_solver_mode
label: _telemetry_solver_mode
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _telemetry_solver_mode
  lines:
  - 46
  - 46
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
  type: String
  units: n/a
  description: Return value of `_telemetry_solver_mode`.
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

# _telemetry_solver_mode

## Purpose
Selects the solver mode string for verification runs from a telemetry-specific variable, then the global one, then a default.

## Design & Implementation
Reads `SPACEAGORA_TELEMETRY_SOLVER_MODE`, falls back to `SPACEAGORA_SOLVER_MODE`, and returns `auto_stiff` if both are empty. `@inline` with a `::String` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_telemetry_solver_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_once|_run_once]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:35-35`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The mode string is not validated here, so a typo reaches the solver policy layer before being rejected.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 46.
