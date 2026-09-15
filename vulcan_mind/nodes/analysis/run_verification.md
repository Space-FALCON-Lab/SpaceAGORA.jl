---
id: analysis.run_verification
label: run_verification
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: run_verification
  lines:
  - 478
  - 481
inputs:
- id: request
  type: VerificationRequest
  units: n/a
  required: true
  description: Positional argument `request`.
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying the request and simulation
    services.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: VerificationResult
  units: n/a
  description: Return value of `run_verification`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
- verification
charts:
- analysis
origin: agent
---

# run_verification

## Purpose
`run_verification` is the public runner for one telemetry verification request. It converts the request into a study configuration, executes the configured scenario through the simulation engine, loads the expected telemetry, computes comparison metrics, and returns a `VerificationResult` for reporting or automated regression checks.

## Theory & Math
Each observed quantity is compared using a combined absolute and relative bound, `|actual-reference| ≤ atol + rtol|reference|`. Orbital and atmospheric channels use separate tolerances and timestep policies defined by the parent module. Diagnostic tables preserve the residuals that explain a failed comparison.

## Model & Assumptions
The request must identify a supported scenario and provide paths or defaults that resolve to valid telemetry inputs. The simulation and reference data must share units, frames, epoch conventions, and sample alignment. Verification thresholds describe acceptance criteria; they do not estimate physical uncertainty.

## Design & Implementation
`runner.jl` defines `run_verification` as a thin public wrapper around `_run_verification(_study_config(request))`. The private runner builds the scenario, invokes `SimulationEngine.run_simulation`, loads reference values, calls comparison and decay-diagnostic helpers, and packages reporting data. CLI and study entrypoints reuse the same result path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `request` | VerificationRequest | n/a | yes | Positional argument `request`. |
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying the request and simulation services. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VerificationResult | n/a | — | Return value of `run_verification`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.ic_fit__ic_fit_run|_ic_fit_run]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:59-59`
- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · feedback · `src/analysis/verification/telemetry_verification/runner.jl:474-474`
- [[analysis.runner_run_verification_cli|run_verification_cli]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:500-500`

**Downstream**

- `callees` → [[analysis.manifest_parsing__study_config|_study_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:479-479`
- `callees` → [[analysis.runner__run_verification|_run_verification]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:479-479`
<!-- vulcan:connections:end -->

## Limitations
Missing telemetry, incompatible SPICE data, solver failures, and scenario-construction errors can prevent a result or mark it failed. Passing a tolerance check only establishes agreement with the selected reference files. A result may be sensitive to solver tolerances, interpolation, telemetry decimation, and the exact manifest revision.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl:478-481`.
