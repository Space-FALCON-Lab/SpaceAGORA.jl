---
id: analysis.runner_run_study
label: run_study
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: run_study
  lines:
  - 509
  - 509
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
  type: Tuple
  units: n/a
  description: Return value of `run_study`. Returns `(summary=result.summary, errors=result.errors)`.
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

# run_study

## Purpose
`run_study` is a zero-argument convenience wrapper for script-style entrypoints that runs the telemetry verification CLI flow on the process `ARGS` and returns only the two result tables. It exists so legacy scripts can `include` the runner and call one function.

## Design & Implementation
Signature `run_study()`. It calls `run_verification_cli(copy(ARGS))` and returns the `NamedTuple` `(summary=result.summary, errors=result.errors)`, discarding the `VerificationResult`'s paths, `plots_dir`, `profile`, `enforce` flag and `total_runtime_s`. Because it delegates entirely, all argument parsing, scenario selection, simulation, CSV writing and threshold enforcement behave exactly as in `run_verification_cli`; a threshold failure under `--enforce` surfaces as the same thrown error. No arguments can be injected programmatically.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `run_study`. Returns `(summary=result.summary, errors=result.errors)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner_run_verification_cli|run_verification_cli]] · `callees` → `callers` · feedback · `src/analysis/verification/telemetry_verification/runner.jl:504-504`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callees` → [[analysis.runner_run_verification_cli|run_verification_cli]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:510-510`
<!-- vulcan:connections:end -->

## Limitations
Reading `ARGS` makes the function untestable without manipulating the global argument vector. The returned tuple drops provenance (output paths, plots directory, runtime), so callers needing those must use `run_verification_cli` or `run_verification` instead. The name suggests a generic study but it is hard-wired to the telemetry verification workflow.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 509.
