---
id: core.simulation_configuration_simulationsettings
label: SimulationSettings
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: SimulationSettings
  lines:
  - 124
  - 124
inputs:
- id: results
  type: Bool
  units: n/a
  required: false
  description: Field `results` (default `true`).
- id: verbose
  type: Bool
  units: n/a
  required: false
  description: Field `verbose` (default `false`).
- id: results_directory
  type: String
  units: n/a
  required: false
  description: Field `results_directory` (default `"output"`).
- id: generate_plots
  type: Bool
  units: n/a
  required: false
  description: Field `generate_plots` (default `true`).
- id: generate_filenames
  type: Bool
  units: n/a
  required: false
  description: Field `generate_filenames` (default `false`).
- id: normalize
  type: Bool
  units: n/a
  required: false
  description: Field `normalize` (default `false`).
- id: save_csv
  type: Bool
  units: n/a
  required: false
  description: Field `save_csv` (default `true`).
- id: checkpoint_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `checkpoint_enabled` (default `false`).
- id: checkpoint_interval_s
  type: Float64
  units: n/a
  required: false
  description: Field `checkpoint_interval_s` (default `300.0`).
- id: checkpoint_directory
  type: String
  units: n/a
  required: false
  description: Field `checkpoint_directory` (default `""`).
- id: resume_from_checkpoint
  type: Bool
  units: n/a
  required: false
  description: Field `resume_from_checkpoint` (default `false`).
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
  type: SimulationSettings
  units: n/a
  description: Constructed `SimulationSettings` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# SimulationSettings

## Purpose
`SimulationSettings` gathers the run-level switches that control output, logging, plotting, and checkpoint/restart behaviour without affecting the physics. It is a field of `SimulationConfiguration` and is rewritten wholesale by the telemetry verification runner to force CSV output into a temporary directory.

## Design & Implementation
A `@kwdef struct` with eleven fields: `results::Bool = true`, `verbose::Bool = false`, `results_directory::String = "output"`, `generate_plots::Bool = true`, `generate_filenames::Bool = false` (whether output names embed run parameters), `normalize::Bool = false` (retained for legacy compatibility; typed runs propagate SI state directly), `save_csv::Bool = true` (CSV in addition to Feather), `checkpoint_enabled::Bool = false`, `checkpoint_interval_s::Float64 = 300.0` seconds of simulated time, `checkpoint_directory::String = ""` (empty means `results_directory/checkpoints`), and `resume_from_checkpoint::Bool = false`. All fields are plain data with no validation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results` | Bool | n/a | no | Field `results` (default `true`). |
| in | `verbose` | Bool | n/a | no | Field `verbose` (default `false`). |
| in | `results_directory` | String | n/a | no | Field `results_directory` (default `"output"`). |
| in | `generate_plots` | Bool | n/a | no | Field `generate_plots` (default `true`). |
| in | `generate_filenames` | Bool | n/a | no | Field `generate_filenames` (default `false`). |
| in | `normalize` | Bool | n/a | no | Field `normalize` (default `false`). |
| in | `save_csv` | Bool | n/a | no | Field `save_csv` (default `true`). |
| in | `checkpoint_enabled` | Bool | n/a | no | Field `checkpoint_enabled` (default `false`). |
| in | `checkpoint_interval_s` | Float64 | n/a | no | Field `checkpoint_interval_s` (default `300.0`). |
| in | `checkpoint_directory` | String | n/a | no | Field `checkpoint_directory` (default `""`). |
| in | `resume_from_checkpoint` | Bool | n/a | no | Field `resume_from_checkpoint` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationSettings | n/a | — | Constructed `SimulationSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:36-36`
- [[analysis.runner__run_simulation_dataframe|_run_simulation_dataframe]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:10-10`
- [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:663-663`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:136-136`
- [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:237-237`
- [[simulation.constellation_ensemble__ensemble_member_settings|_ensemble_member_settings]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:15-15`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`checkpoint_interval_s` is not checked for positivity and a zero value would checkpoint every step. `normalize` is a no-op flag whose presence can mislead readers into expecting state scaling. `results_directory` duplicates the intent of `FilePaths.results`. Setting `results=false` with `save_csv=true` is not rejected, so the meaning of that combination is defined only by the engine.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 124.
