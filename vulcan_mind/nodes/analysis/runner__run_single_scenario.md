---
id: analysis.runner__run_single_scenario
label: _run_single_scenario
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: _run_single_scenario
  lines:
  - 134
  - 134
inputs:
- id: cfg
  type: OrbitEventsScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: Any
  units: n/a
  description: Return value of `_run_single_scenario`. Returns `annotated_rows, final_errors,
    calibration_runtime_s`.
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

# _run_single_scenario

## Purpose
`_run_single_scenario` runs one manifest scenario end to end (optional calibration grid search followed by a final evaluation) and returns annotated metric rows, per-event error tables, and the scenario's runtime. Two methods exist, for `OrbitEventsScenarioConfig` (apoapsis/periapsis event comparison) and `TimeAlignedScenarioConfig` (sample-by-sample comparison against telemetry).

## Design & Implementation
Both take `(cfg, profile::Symbol)`. They compute `use_calibration = _calibration_active(cfg.calibration, profile)`, choose quick/full orbit counts or point caps, and, when `cal.search_on_quick_subset` is set with a `:full` profile, evaluate the grid on the `:quick` subset. Candidate grids come from `_grid_values(min, max, steps)` for `cd_scale` (if `fit_cd_scale`) and `cr` (if `fit_cr && srp_enabled`), else the single defaults `[1.0]` and `[cfg.srp_cr]`. For each `(cd_scale, cr)` pair they build args via `_make_orbit_args` or `_make_time_aligned_args`, apply `_with_study_settings`, run `_run_simulation_dataframe`, compute rows and errors (`_orbit_rows_errors` or `_time_aligned_rows_errors`), optionally re-evaluate with `_estimate_event_biases`, and keep the minimum `_calibration_score`. The final run uses the best pair and `_final_run_or_reused_eval` to avoid a redundant solve; results are annotated by `_annotate_calibration_rows` with `best_cd`, `best_cr`, biases, score, runtimes, `dt_max_orbit` and solver info. The time-aligned method additionally loads telemetry, derives the initial condition, and sets mission time to `max(time_s[end], 1.0)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_run_single_scenario`. Returns `annotated_rows, final_errors, calibration_runtime_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:378-378`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callees` → [[analysis.calibration__annotate_calibration_rows|_annotate_calibration_rows]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:197-197`
- `callees` → [[analysis.calibration__calibration_active|_calibration_active]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:137-137`
- `callees` → [[analysis.calibration__calibration_score|_calibration_score]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:171-171`
- `callees` → [[analysis.calibration__grid_values|_grid_values]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:149-149`
- `callees` → [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:257-257`
- `callees` → [[analysis.runner__final_run_or_reused_eval|_final_run_or_reused_eval]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:184-184`
- `callees` → [[analysis.runner__initial_condition_from_time_aligned_telemetry|_initial_condition_from_time_aligned_telemetry]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:220-220`
- `callees` → [[analysis.runner__run_simulation_dataframe|_run_simulation_dataframe]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:163-163`
- `callees` → [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:161-161`
- `callees` → [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:246-246`
- `callees` → [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:162-162`
- `callees` → [[analysis.telemetry_loading__load_time_aligned_telemetry|_load_time_aligned_telemetry]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:219-219`
- `callees` → [[envana.ana_calibration_estimate_event_biases|_estimate_event_biases]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:168-168`
- `callees` → [[envana.ana_error_tables_orbit_rows_errors|_orbit_rows_errors]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:166-166`
<!-- vulcan:connections:end -->

## Limitations
The grid search is exhaustive with cost `cd_steps × cr_steps` full simulations and no early termination. `best_score` starts at `Inf`, so if every score is `NaN` the defaults are kept silently. Ties keep the first candidate in iteration order. The two methods duplicate roughly 60 lines of identical control flow, so fixes must be applied twice. `reused_eval_run` always holds the last grid point evaluated, not the best one; correctness relies on `_single_point_calibration` only permitting reuse when the grid has one point.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 134.
