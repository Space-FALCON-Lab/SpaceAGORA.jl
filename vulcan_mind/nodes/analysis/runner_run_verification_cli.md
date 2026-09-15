---
id: analysis.runner_run_verification_cli
label: run_verification_cli
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: run_verification_cli
  lines:
  - 488
  - 488
inputs:
- id: args
  type: Vector{String}
  units: n/a
  required: false
  description: Positional argument `args` (default `copy(ARGS)`).
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
  type: VerificationResult
  units: n/a
  description: Return value of `run_verification_cli`.
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

# run_verification_cli

## Purpose
`run_verification_cli` is the command-line entrypoint for the telemetry orbit-accuracy study: it parses `ARGS`-style options into a `StudyConfig`, converts that into a typed `VerificationRequest`, and executes the study. It is exported from the package and invoked by the CLI `telemetry` command and by direct script execution.

## Design & Implementation
Signature `run_verification_cli(args::Vector{String}=copy(ARGS))::VerificationResult`. It calls `parse_cli(args)` to obtain `cfg`, then `_request_from_study_config(cfg)` for a base request, and rebuilds a `VerificationRequest` copying `profile`, `out_summary`, `out_errors` and `manifest_path` from that base while taking `enforce`, `generate_plots` and `scenarios` directly from `cfg`. The rebuilt request is passed to `run_verification`, which in turn calls `_run_verification(_study_config(request))`. A trailing `if abspath(PROGRAM_FILE) == abspath(@__FILE__)` guard at file scope runs this function when the file is executed as a script.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vector{String} | n/a | no | Positional argument `args` (default `copy(ARGS)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VerificationResult | n/a | — | Return value of `run_verification_cli`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner_run_study|run_study]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:510-510`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__request_from_study_config|_request_from_study_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:490-490`
- `callees` → [[analysis.run_verification|run_verification]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:500-500`
- `callees` → [[analysis.runner_run_study|run_study]] · `callers` · feedback · `src/analysis/verification/telemetry_verification/runner.jl:504-504`
- `callees` → [[analysis.types_verificationrequest|VerificationRequest]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:491-491`
- `callees` → [[envana.ana_manifest_parsing_parse_cli|parse_cli]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:489-489`
<!-- vulcan:connections:end -->

## Limitations
The request is round-tripped `StudyConfig -> VerificationRequest -> StudyConfig`, so any field that `_request_from_study_config` does not carry is silently dropped unless explicitly re-copied here (only three fields are). Parse errors from `parse_cli` propagate as exceptions rather than usage messages. The default `copy(ARGS)` captures the Julia process arguments, which in a REPL or test harness may contain unrelated flags.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 488.
