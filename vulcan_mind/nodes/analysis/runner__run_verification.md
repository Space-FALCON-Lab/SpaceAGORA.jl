---
id: analysis.runner__run_verification
label: _run_verification
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: _run_verification
  lines:
  - 354
  - 354
inputs:
- id: cfg
  type: StudyConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_run_verification`.
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

# _run_verification

## Purpose
`_run_verification` is the study driver: it loads and selects scenarios, runs each one, merges the metric rows with threshold gates and provenance columns, writes the summary and error CSVs, optionally generates plots, and enforces pass/fail thresholds. It returns a `VerificationResult`.

## Design & Implementation
Signature `(cfg::StudyConfig)::VerificationResult`. It calls `_select_scenarios(_load_scenarios_from_manifest(cfg.manifest_path), cfg.scenarios)` and prints a header (profile, enforce flag, manifest path, deterministic GRAM mode). For each scenario it prints the atmosphere-truth settings, calls `_run_single_scenario(sc, cfg.profile)`, then for each returned row merges it with `_evaluate_thresholds(row, sc, profile)` and a `NamedTuple` of metadata (source file, units, maneuver counts, `timestamp_utc`, all `atmosphere_truth` fields including `gram_seed` and `gram_perturbation_scales`). Error tables get a `dt_max_orbit_s` column looked up by event. After the loop it builds `summary_df = DataFrame(summary_rows)`, `errors_df` via `vcat(...; cols=:union)`, sorts by `[:scenario, :event]`, adds `total_runtime_s`, calls `_append_display_metric_columns!` and `_append_display_error_columns!`, creates output directories with `mkpath`, and writes both CSVs. Plots are produced by `_generate_plots` when `cfg.generate_plots`. When `cfg.enforce` and any `pass == false`, it shows the failing rows and calls `error`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | StudyConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VerificationResult | n/a | — | Return value of `_run_verification`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.run_verification|run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:479-479`

**Downstream**

- `callees` → [[analysis.calibration__calibration_active|_calibration_active]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:376-376`
- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:412-412`
- `callees` → [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:355-355`
- `callees` → [[analysis.reporting__append_display_error_columns_bang|_append_display_error_columns!]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:432-432`
- `callees` → [[analysis.reporting__append_display_metric_columns_bang|_append_display_metric_columns!]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:431-431`
- `callees` → [[analysis.reporting__axis_units|_axis_units]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:389-389`
- `callees` → [[analysis.reporting__generate_plots|_generate_plots]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:446-446`
- `callees` → [[analysis.reporting__maneuver_count|_maneuver_count]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:392-392`
- `callees` → [[analysis.reporting__maneuver_replay_scale_mode|_maneuver_replay_scale_mode]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:393-393`
- `callees` → [[analysis.reporting__orbit_altitude_mode|_orbit_altitude_mode]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:391-391`
- `callees` → [[analysis.reporting__scenario_status_extra|_scenario_status_extra]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:376-376`
- `callees` → [[analysis.reporting__source_file|_source_file]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:388-388`
- `callees` → [[analysis.reporting__value_units|_value_units]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:390-390`
- `callees` → [[analysis.run_verification|run_verification]] · `callers` · feedback · `src/analysis/verification/telemetry_verification/runner.jl:474-474`
- `callees` → [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:378-378`
- `callees` → [[analysis.runner__select_scenarios|_select_scenarios]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:355-355`
- `callees` → [[analysis.types_verificationresult|VerificationResult]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:461-461`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:361-361`
- `callees` → [[envana.ana_reporting_evaluate_thresholds|_evaluate_thresholds]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:382-382`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:383-383`
<!-- vulcan:connections:end -->

## Limitations
Scenarios run strictly sequentially; there is no parallelism across scenarios. `DataFrame(summary_rows)` requires every merged row to have identical keys, so a scenario type that yields different metric columns causes a construction error at the end after all simulations have run. Threshold failure is reported by throwing after files are already written, so callers must inspect the CSVs for detail. Console output uses `println`/`show` directly and cannot be redirected through a logger. An empty scenario list produces empty frames and skips sorting without warning.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 354.
