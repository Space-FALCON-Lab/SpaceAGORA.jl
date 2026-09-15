---
id: analysis.runner__run_once
label: _run_once
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: _run_once
  lines:
  - 32
  - 32
inputs:
- id: maxiters
  type: Int
  units: n/a
  required: true
  description: Positional argument `maxiters`.
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
  type: Any
  units: n/a
  description: Return value of `_run_once`. Returns `solve_result, elapsed_s`.
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

# _run_once

## Purpose
`_run_once` is the closure inside `_run_simulation_dataframe` that performs a single simulation attempt with a given solver iteration cap. It exists so the enclosing function can call the same fully-configured run twice: once with the profile's default `maxiters` and, on a MaxIters failure, once with the enlarged retry value.

## Design & Implementation
Signature `_run_once(maxiters::Int)`, capturing `tmp`, `cfg_run`, `save_fields` and `truth` from the enclosing scope. It measures the run with `@elapsed`, reads `solver_mode = _telemetry_solver_mode()`, and wraps the call in `withenv` setting `SPACEAGORA_WARN_NORMALIZE=0`, `SPACEAGORA_WARN_DEPRECATED_CONFIG=0`, `SPACEAGORA_SOLVER_MODE`, `SPACEAGORA_SOLVER_MAXITERS=string(maxiters)`, `SPACEAGORA_GRAM_OFFLINE_SURROGATE`, `SPACEAGORA_GRAM_STATIC_GRID` and `SPACEAGORA_GRAM_TRACK_CACHE` (`on`/`off` from the truth booleans) and `SPACEAGORA_GRAM_GLOBAL_LOCK`. Inside `cd(tmp)` it calls `SimulationEngine.run_simulation(cfg_run; isolate_state=false, save_fields, return_solution=true, return_solver_metadata=true)` and returns `(solve_result, elapsed_s)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `maxiters` | Int | n/a | yes | Positional argument `maxiters`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_run_once`. Returns `solve_result, elapsed_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__telemetry_solver_mode|_telemetry_solver_mode]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:35-35`
- `callees` → [[analysis.manifest_parsing__telemetry_solver_retry_maxiters|_telemetry_solver_retry_maxiters]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:69-69`
- `callees` → [[simulation.run_simulation|run_simulation]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:47-47`
- `callees` → [[simx.engine_execution_run_simulation|run_simulation]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:47-47`
- `callees` → [[spaceagora.run_simulation|run_simulation]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
It changes the process working directory and environment for the duration of the call; both are process-global, so running two scenarios on different threads at once would race. `withenv` restores variables afterwards but does not protect against the engine caching environment-derived settings on first use. `solve_result` stays `nothing` if the engine throws, and the closure does not catch anything itself. Passing `isolate_state=false` means a mutated `cfg_run` from the first attempt is reused by the retry.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 32.
