---
id: analysis.calibration__single_point_calibration
label: _single_point_calibration
kind: function
source:
  file: src/analysis/verification/telemetry_verification/calibration.jl
  symbol: _single_point_calibration
  lines:
  - 13
  - 13
inputs:
- id: use_calibration
  type: Bool
  units: n/a
  required: true
  description: Positional argument `use_calibration`.
- id: cd_candidates
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `cd_candidates`.
- id: cr_candidates
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `cr_candidates`.
- id: eval_profile
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `eval_profile`.
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
  type: Bool
  units: n/a
  description: Return value of `_single_point_calibration`.
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

# _single_point_calibration

## Purpose

`_single_point_calibration` detects the degenerate calibration case where the search grid collapses to a single point and the evaluation solve therefore has exactly the same configuration as the final solve. When it returns `true` the runner reuses the evaluation result instead of solving the trajectory a second time.

## Design & Implementation

The `@inline`d predicate returns `use_calibration && length(cd_candidates) == 1 && length(cr_candidates) == 1 && eval_profile == profile`. All four conditions must hold: calibration must be active, both the drag-scale and reflectivity candidate vectors must contain exactly one entry, and the evaluation profile symbol must equal the final profile symbol. The docstring records the important consequence: the bias fit is still applied to the final rows either way, so skipping the duplicate solve changes runtime only, never the reported numbers.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `use_calibration` | Bool | n/a | yes | Positional argument `use_calibration`. |
| in | `cd_candidates` | AbstractVector | n/a | yes | Positional argument `cd_candidates`. |
| in | `cr_candidates` | AbstractVector | n/a | yes | Positional argument `cr_candidates`. |
| in | `eval_profile` | Symbol | n/a | yes | Positional argument `eval_profile`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_single_point_calibration`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.calibration__calibration_active|_calibration_active]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:6-6`
- [[analysis.runner__final_run_or_reused_eval|_final_run_or_reused_eval]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:327-327`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Equality of profiles is compared by `Symbol` identity, so two profiles that are numerically identical but differently named still force a redundant solve. The function inspects only candidate-vector lengths, not their values, and takes `AbstractVector` so any lazily generated grid is measured by `length` alone. It says nothing about whether the single candidate is a sensible one.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/calibration.jl` line 13.
