---
id: analysis.runner__final_run_or_reused_eval
label: _final_run_or_reused_eval
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: _final_run_or_reused_eval
  lines:
  - 315
  - 315
inputs:
- id: solve
  type: F
  units: n/a
  required: true
  description: Positional argument `solve`.
- id: reused_eval_run
  type: Any
  units: n/a
  required: true
  description: Positional argument `reused_eval_run`.
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
- id: scenario_name
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `scenario_name`.
- id: best_cd
  type: Float64
  units: n/a
  required: true
  description: Positional argument `best_cd`.
- id: best_cr
  type: Float64
  units: n/a
  required: true
  description: Positional argument `best_cr`.
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
  description: 'Return value of `_final_run_or_reused_eval`. Returns `reused_eval_run`
    or `solve()`. Type parameters: `{F <: Function}`.'
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

# _final_run_or_reused_eval

## Purpose
`_final_run_or_reused_eval` decides whether the final scenario solve can reuse the calibration-grid evaluation already performed, avoiding a second identical simulation when the grid collapsed to a single point at the same profile. It is invoked with a `do`-block by both `_run_single_scenario` methods.

## Design & Implementation
Signature `(solve::F, reused_eval_run, use_calibration::Bool, cd_candidates::AbstractVector, cr_candidates::AbstractVector, eval_profile::Symbol, profile::Symbol, scenario_name::AbstractString, best_cd::Float64, best_cr::Float64) where {F <: Function}`. If `reused_eval_run !== nothing` and `_single_point_calibration(use_calibration, cd_candidates, cr_candidates, eval_profile, profile)` returns true, it prints a message naming the scenario, `cd_scale` and `cr`, and returns `reused_eval_run` unchanged. Otherwise it calls the zero-argument `solve()` closure and returns its result. It has no side effects beyond the `println`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `solve` | F | n/a | yes | Positional argument `solve`. |
| in | `reused_eval_run` | Any | n/a | yes | Positional argument `reused_eval_run`. |
| in | `use_calibration` | Bool | n/a | yes | Positional argument `use_calibration`. |
| in | `cd_candidates` | AbstractVector | n/a | yes | Positional argument `cd_candidates`. |
| in | `cr_candidates` | AbstractVector | n/a | yes | Positional argument `cr_candidates`. |
| in | `eval_profile` | Symbol | n/a | yes | Positional argument `eval_profile`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `scenario_name` | AbstractString | n/a | yes | Positional argument `scenario_name`. |
| in | `best_cd` | Float64 | n/a | yes | Positional argument `best_cd`. |
| in | `best_cr` | Float64 | n/a | yes | Positional argument `best_cr`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_final_run_or_reused_eval`. Returns `reused_eval_run` or `solve()`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:184-184`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callees` → [[analysis.calibration__single_point_calibration|_single_point_calibration]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:327-327`
- `callees` → [[analysis.runner__select_scenarios|_select_scenarios]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:335-335`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:328-328`
<!-- vulcan:connections:end -->

## Limitations
Reuse is justified only because the caller guarantees the single grid point's configuration equals the final one; the function does not verify that `args_eval` and `args_final` match (for example that orbit counts or point caps agree), so a caller change could make reuse incorrect. `reused_eval_run` is untyped and returned as-is. Progress reporting goes to `stdout` via `println` rather than a logger, so it cannot be silenced.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 315.
