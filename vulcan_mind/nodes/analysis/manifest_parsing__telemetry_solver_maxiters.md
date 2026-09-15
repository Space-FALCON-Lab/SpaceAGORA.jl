---
id: analysis.manifest_parsing__telemetry_solver_maxiters
label: _telemetry_solver_maxiters
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _telemetry_solver_maxiters
  lines:
  - 35
  - 35
inputs:
- id: profile
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `profile`.
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
  description: Return value of `_telemetry_solver_maxiters`.
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

# _telemetry_solver_maxiters

## Purpose
Resolves the integrator iteration cap for a verification profile, allowing an environment override.

## Design & Implementation
Chooses `TELEMETRY_SOLVER_MAXITERS_QUICK_DEFAULT` or `_FULL_DEFAULT` by profile, then reads `SPACEAGORA_TELEMETRY_SOLVER_MAXITERS` through `_parse_positive_int_env` with that default. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_telemetry_solver_maxiters`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_simulation_dataframe|_run_simulation_dataframe]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:30-30`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_int_env|_parse_positive_int_env]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:37-37`
<!-- vulcan:connections:end -->

## Limitations
One environment variable overrides both profiles, so a quick-tuned override applied in a full run can be far too small.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 35.
