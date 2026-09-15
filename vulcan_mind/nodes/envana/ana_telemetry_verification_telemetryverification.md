---
id: envana.ana_telemetry_verification_telemetryverification
label: TelemetryVerification
kind: struct
source:
  file: src/analysis/verification/telemetry_verification.jl
  symbol: TelemetryVerification
  lines:
  - 1
  - 48
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Analysis namespace root that declares the verification exports and
    includes its twelve implementation files.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: verification_api
  type: Module
  units: n/a
  description: Exported VerificationRequest, VerificationResult, run_verification,
    run_verification_cli and run_study surface.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- envana
origin: agent
---
# TelemetryVerification

## Purpose
`TelemetryVerification` is the module that compares simulated trajectories against recorded flight or reference telemetry and decides whether a scenario passes. It declares the strict-profile numerical constants, the default manifest and output paths, and then includes the twelve files that implement parsing, scenario construction, telemetry loading, metric computation, calibration, error tabulation, reporting, execution, and initial-condition fitting.

## Theory & Math
The strict profile fixes the integrator error budget used when a run must be reproducible to reference accuracy: relative and absolute orbit tolerances of `1e-7` and `1e-9` (dimensionless and in state units respectively), the same pair for the atmospheric phase, a maximum orbit step of 60.0 s, and a maximum atmospheric step of 0.2 s. The step ceiling matters because local truncation error for a Runge-Kutta method of order `p` scales as `O(dt^(p+1))`, so capping `dt` bounds per-step error independently of the adaptive controller. Solver iteration ceilings of 5,000,000 for the quick profile and 20,000,000 for the full profile bound wall-clock cost.

## Model & Assumptions
`REPO_ROOT` is derived by walking three directories up from `@__DIR__`, so the module assumes it sits at `src/analysis/verification/`. Outputs default to `output/` and the manifest defaults to `test/telemetry_benchmark_manifest.toml`. SPICE kernels are expected under the GRAM Suite SPICE directory named by `SPICE_PATH`. The trailing comment states the architecture contract explicitly: the runner delegates execution to `SimulationEngine.run_simulation` rather than driving the integrator itself.

## Design & Implementation
Include order is deliberate and matches the dependency chain: `types.jl` first so configuration structs exist, then `manifest_parsing.jl`, `example_support.jl`, `scenario_builders.jl`, `telemetry_loading.jl`, `comparison_metrics.jl`, `decay_diagnostics.jl`, `calibration.jl`, `error_tables.jl`, `reporting.jl`, `runner.jl`, and `ic_fit.jl`. Data handling relies on `CSV`, `DataFrames`, and `Arrow`; configuration on `TOML`; numerics on `Statistics`, `LinearAlgebra`, and `StaticArrays`. A `reference_system.jl` include from `src/core/interfaces` supplies frame definitions.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Analysis namespace root that declares the verification exports and includes its twelve implementation files. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `verification_api` | Module | n/a | — | Exported VerificationRequest, VerificationResult, run_verification, run_verification_cli and run_study surface. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because paths are computed relative to the source file, relocating the module breaks manifest and kernel discovery unless the environment variables are set. Only five names are exported, so every internal metric and builder is private and cannot be reused by external analysis code without qualification. The strict constants are compile-time `const` values and cannot be tuned per scenario.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification.jl:1-48`, covering the constant block, the dependency imports, and the full include list.
