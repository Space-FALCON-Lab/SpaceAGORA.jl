---
id: analysis.manifest_parsing__telemetry_solver_retry_maxiters
label: _telemetry_solver_retry_maxiters
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _telemetry_solver_retry_maxiters
  lines:
  - 40
  - 40
inputs:
- id: base_maxiters
  type: Int
  units: n/a
  required: true
  description: Positional argument `base_maxiters`.
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
  description: Return value of `_telemetry_solver_retry_maxiters`.
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

# _telemetry_solver_retry_maxiters

## Purpose
Chooses the iteration cap for the automatic retry after a MaxIters failure, guaranteeing it exceeds the first attempt.

## Design & Implementation
Defaults to the larger of four times the base cap and the base plus one million, reads `SPACEAGORA_TELEMETRY_SOLVER_MAXITERS_RETRY` with that default, and finally clamps the result to at least `base_maxiters + 1`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `base_maxiters` | Int | n/a | yes | Positional argument `base_maxiters`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_telemetry_solver_retry_maxiters`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_once|_run_once]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:69-69`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_int_env|_parse_positive_int_env]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:42-42`
<!-- vulcan:connections:end -->

## Limitations
The retry is a single fixed enlargement; a run that needs more than that fails outright rather than escalating again.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 40.
