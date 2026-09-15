---
id: analysis.example_support__example_smoke_args
label: _example_smoke_args
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: _example_smoke_args
  lines:
  - 19
  - 19
inputs:
- id: args
  type: SM.SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: SM.SimulationConfiguration
  units: n/a
  description: Return value of `_example_smoke_args`. Returns `args` or `SM.SimulationConfiguration(`.
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

# _example_smoke_args

## Purpose
Rewrites a full example configuration into its smoke-test equivalent — one orbit, a short mission, few saved steps, plots off and results local — leaving the physics models untouched.

## Design & Implementation
Returns `args` unchanged unless `_example_smoke_enabled()`. Otherwise it constructs a new `MissionConfiguration` keeping `mission_type`, `keplerian`, `orientation_sim` and `data_rate` but forcing `number_of_orbits` to 1, `mission_time` through `_example_smoke_mission_time`, and `num_steps_to_save` clamped into the range 50 to 200. A new `SimulationSettings` disables verbosity, plotting and normalisation, ties `results` and `save_csv` to the results-enabled flag, and points `results_directory` at `pwd()/output` so concurrent smoke runs in different directories cannot collide. Every other field of the `SimulationConfiguration` is passed through by name.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SM.SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SM.SimulationConfiguration | n/a | — | Return value of `_example_smoke_args`. Returns `args` or `SM.SimulationConfiguration(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support_run_and_report|run_and_report]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:178-178`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`

**Downstream**

- `callees` → [[analysis.example_support__example_smoke_enabled|_example_smoke_enabled]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:20-20`
- `callees` → [[analysis.example_support__example_smoke_mission_time|_example_smoke_mission_time]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:31-31`
- `callees` → [[analysis.example_support__example_smoke_results_enabled|_example_smoke_results_enabled]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:26-26`
- `callees` → [[core.simulation_configuration_missionconfiguration|MissionConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:27-27`
- `callees` → [[core.simulation_configuration_simulationsettings|SimulationSettings]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:36-36`
- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
The pass-through is by explicit field list, so a field added to `SimulationConfiguration` later is dropped in smoke mode until this function is updated; `solver_config` in particular is not forwarded, so a smoke run always uses the default solver policy.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 19.
