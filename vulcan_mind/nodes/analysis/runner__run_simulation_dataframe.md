---
id: analysis.runner__run_simulation_dataframe
label: _run_simulation_dataframe
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: _run_simulation_dataframe
  lines:
  - 1
  - 1
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: scenario_name
  type: String
  units: n/a
  required: true
  description: Positional argument `scenario_name`.
- id: truth
  type: AtmosphereTruthConfig
  units: n/a
  required: true
  description: Positional argument `truth`.
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
  type: Tuple
  units: n/a
  description: Return value of `_run_simulation_dataframe`. Returns `mktempdir() do
    tmp` or `solve_result, elapsed_s` or `(results_df=results_df, elapsed_s=elapsed_s,
    solver_info=solver_info)`.
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

# _run_simulation_dataframe

## Purpose
`_run_simulation_dataframe` executes one telemetry-verification simulation in an isolated temporary directory and returns its CSV results as a `DataFrame` together with wall-clock time and solver metadata. Both `_run_single_scenario` methods call it for every calibration grid point and for the final run.

## Design & Implementation
Signature `(args::SimulationConfiguration, scenario_name::String, truth::AtmosphereTruthConfig, profile::Symbol)`. Inside `mktempdir() do tmp` it rebuilds `cfg_run` as a copy of `args` whose `SimulationSettings` force `results=true`, `verbose=false`, `results_directory=tmp`, `generate_plots=false`, `generate_filenames=false`, `normalize=false`, `save_csv=true`. It fetches `save_fields = _save_fields_for_study()` and `base_maxiters = _telemetry_solver_maxiters(profile)`, then defines the closure `_run_once(maxiters)` and calls it; if the raised error satisfies `_is_maxiters_error` it retries once with `_telemetry_solver_retry_maxiters(base_maxiters)` and logs a `@warn`, otherwise rethrows. After the run it requires `tmp/simulation_results.csv` to exist (else `error`), reads it with `CSV.read`, and condenses `solve_result.solver_trace` into a `solver_info` `NamedTuple` (`solver_mode`, `solver_sequence` joined by `->`, fallback count and trigger retcodes joined by `|`, `solver_retcode`, `solver_maxiters`, `solver_maxiters_retry_used`). Returns `(results_df, elapsed_s, solver_info)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `scenario_name` | String | n/a | yes | Positional argument `scenario_name`. |
| in | `truth` | AtmosphereTruthConfig | n/a | yes | Positional argument `truth`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_run_simulation_dataframe`. Returns `mktempdir() do tmp` or `solve_result, elapsed_s` or `(results_df=results_df, elapsed_s=elapsed_s, solver_info=solver_info)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:163-163`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__telemetry_solver_maxiters|_telemetry_solver_maxiters]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:30-30`
- `callees` → [[analysis.scenario_builders__save_fields_for_study|_save_fields_for_study]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:29-29`
- `callees` → [[core.simulation_configuration_simulationsettings|SimulationSettings]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:10-10`
- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:8-8`
<!-- vulcan:connections:end -->

## Limitations
The temporary directory is deleted when the `do` block exits, so raw simulation artifacts other than the returned frame are lost. Only one retry on MaxIters is attempted; a second failure propagates. The CSV filename `simulation_results.csv` is hard-coded and must match what `IOOutputs` writes with `generate_filenames=false`. Reading results back through CSV loses column type information (for example `Symbol` or `Bool` columns become strings). `isolate_state=false` is passed, so `cfg_run` must not be shared concurrently.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 1.
