---
id: analysis.types_verificationresult
label: VerificationResult
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: VerificationResult
  lines:
  - 276
  - 276
inputs:
- id: summary
  type: DataFrame
  units: n/a
  required: true
  description: Field `summary`.
- id: errors
  type: DataFrame
  units: n/a
  required: true
  description: Field `errors`.
- id: summary_path
  type: String
  units: n/a
  required: true
  description: Field `summary_path`.
- id: errors_path
  type: String
  units: n/a
  required: true
  description: Field `errors_path`.
- id: plots_dir
  type: String
  units: n/a
  required: true
  description: Field `plots_dir`.
- id: profile
  type: Symbol
  units: n/a
  required: true
  description: Field `profile`.
- id: enforce
  type: Bool
  units: n/a
  required: true
  description: Field `enforce`.
- id: total_runtime_s
  type: Float64
  units: n/a
  required: true
  description: Field `total_runtime_s`.
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
  description: Constructed `VerificationResult` (keyword constructor via @kwdef).
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

# VerificationResult

## Purpose
Structured return value of a telemetry verification run, bundling the per-scenario summary table, the pointwise error table, the paths where those tables and plots were written, the profile and enforcement flags used, and total wall-clock runtime.

## Design & Implementation
`Base.@kwdef struct` with no defaults: `summary::DataFrame` and `errors::DataFrame` hold the in-memory tables that were also written to `summary_path` and `errors_path`; `plots_dir::String` points at the plot output directory; `profile::Symbol` and `enforce::Bool` echo the request so callers (for example the nightly harness) can interpret pass/fail columns; `total_runtime_s::Float64` is seconds. Fields are immutable but the `DataFrame` contents are mutable references.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `summary` | DataFrame | n/a | yes | Field `summary`. |
| in | `errors` | DataFrame | n/a | yes | Field `errors`. |
| in | `summary_path` | String | n/a | yes | Field `summary_path`. |
| in | `errors_path` | String | n/a | yes | Field `errors_path`. |
| in | `plots_dir` | String | n/a | yes | Field `plots_dir`. |
| in | `profile` | Symbol | n/a | yes | Field `profile`. |
| in | `enforce` | Bool | n/a | yes | Field `enforce`. |
| in | `total_runtime_s` | Float64 | n/a | yes | Field `total_runtime_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VerificationResult | n/a | — | Constructed `VerificationResult` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:461-461`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct does not carry an explicit pass/fail verdict; callers must inspect the summary columns. Holding both DataFrames in memory duplicates what is on disk and can be large for full-profile runs with many evaluation points. `plots_dir` is populated even when plot generation was disabled, so its existence on disk is not guaranteed. There is no schema check that `summary` and `errors` contain the expected columns.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 276.
