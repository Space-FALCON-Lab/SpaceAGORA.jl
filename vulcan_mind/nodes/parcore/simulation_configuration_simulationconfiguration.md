---
id: parcore.simulation_configuration_simulationconfiguration
label: SimulationConfiguration
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: SimulationConfiguration
  lines:
  - 235
  - 247
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Configuration sub-records (file paths, settings, mission configuration,
    integration tolerances) declared earlier in the same module.
- id: type_params
  type: AbstractType
  units: n/a
  required: true
  description: Abstract planet, density, ephemerides and thermal supertypes that parameterise
    the environment model field.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: config
  type: SimulationConfiguration
  units: n/a
  description: Fully typed, immutable simulation configuration passed to the ODE parameter
    record and to every model evaluation.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- parcore
origin: agent
---

# SimulationConfiguration

## Purpose
`SimulationConfiguration` is the root configuration record of a run. It gathers file paths, general simulation settings, mission configuration, the environment model, the dynamics, guidance, navigation and control model sets, the initial epoch, integrator tolerances and an optional solver configuration into one immutable value that is threaded through the entire propagation.

## Model & Assumptions
The struct is parameterised on five type variables: the planet, density, ephemerides and thermal types are constrained to their abstract supertypes, and the dynamics model set is constrained to a `Tuple`. Making these type parameters rather than abstract fields is what allows the integrator right-hand side to specialise: model dispatch is resolved at compile time instead of through a dynamic lookup on every derivative evaluation.

## Design & Implementation
The record is built with `@kwdef`, and most fields carry defaults, so a caller supplies only the environment and model sets. `solver_config` defaults to `nothing`, which is the explicit signal to read solver choices from the environment at run time rather than from the configuration. Sibling records in the same `SimConfig` module supply the parts: `FilePaths`, `SimulationSettings`, `MissionConfiguration`, `EnvironmentModel`, `InitialTime`, `IntegrationTolerances` and `SolverConfig`. `MissionConfiguration` and `EnvironmentModel` both declare inner constructors that normalise and validate their arguments before the outer record is assembled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Configuration sub-records (file paths, settings, mission configuration, integration tolerances) declared earlier in the same module. |
| in | `type_params` | AbstractType | n/a | yes | Abstract planet, density, ephemerides and thermal supertypes that parameterise the environment model field. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `config` | SimulationConfiguration | n/a | — | Fully typed, immutable simulation configuration passed to the ODE parameter record and to every model evaluation. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:47-47`
- [[analysis.runner__run_simulation_dataframe|_run_simulation_dataframe]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:8-8`
- [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:433-433`
- [[analysis.scenario_builders__with_environment_wind|_with_environment_wind]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:380-380`
- [[analysis.scenario_builders__with_orbit_mission|_with_orbit_mission]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:468-468`
- [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:661-661`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:135-135`
- [[simulation.constellation_ensemble__ensemble_member_configuration|_ensemble_member_configuration]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:37-37`

**Downstream**

- `callees` → [[core.simulation_configuration_filepaths|FilePaths]] · `callers` · call · `src/core/state/simulation_configuration.jl:236-236`
- `callees` → [[core.simulation_configuration_integrationtolerances|IntegrationTolerances]] · `callers` · call · `src/core/state/simulation_configuration.jl:245-245`
- `callees` → [[core.simulation_configuration_missionconfiguration|MissionConfiguration]] · `callers` · call · `src/core/state/simulation_configuration.jl:238-238`
- `callees` → [[core.simulation_configuration_simulationsettings|SimulationSettings]] · `callers` · call · `src/core/state/simulation_configuration.jl:237-237`
<!-- vulcan:connections:end -->

## Limitations
Because the model sets are type parameters, changing a model family produces a different concrete configuration type and therefore a fresh round of compilation, which is a visible cost when sweeping model choices. There is no schema versioning on the record, so a configuration serialised by one revision of the package may not reconstruct against another. Defaults hide omissions: a caller that forgets a mission field silently receives the default rather than an error.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl:235-247`.
