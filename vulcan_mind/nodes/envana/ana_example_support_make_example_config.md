---
id: envana.ana_example_support_make_example_config
label: make_example_config
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: make_example_config
  lines:
  - 119
  - 175
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace re-exporting the SimulationModel configuration
    constructors.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: configuration
  type: SimulationConfiguration
  units: n/a
  description: Fully populated single-vehicle simulation configuration for examples
    and smoke runs.
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
# make_example_config

## Purpose
`make_example_config` assembles a complete `SimulationConfiguration` from a small set of keyword arguments, so that examples and smoke tests can specify only a planet, a spacecraft, a mission duration, and an epoch while still receiving a fully valid configuration tree.

## Theory & Math
The defaults encode a specific numerical contract. Integration tolerances are set to relative and absolute values of `1e-8` for both the orbit and atmospheric phases, with maximum steps of 20.0 s in orbit and 0.2 s in the atmosphere. For an adaptive Runge-Kutta method the controller accepts a step when the estimated local error satisfies `err <= atol + rtol * |y|`, so with `atol = rtol = 1e-8` the absolute floor dominates only for state components far below unity. The entry interface altitude defaults to `EI_km = 300.0` kilometres, and the thermal model uses a Maxwellian accommodation factor of 1.0, meaning fully accommodated diffuse reflection.

## Model & Assumptions
The default dynamics stack is a single `InverseSquaredJ2GravityModel`, the default atmosphere is `NoAtmosphereModel`, and the default ephemerides source is `SpiceEphemeridesModel`, so an example runs as a J2-only vacuum propagation unless overridden. Mission configuration is fixed to `MissionTime` with `number_of_orbits = 1`, `num_steps_to_save = 1000`, and Keplerian initial conditions enabled by default. Guidance, navigation, and control models are constructed empty, with no effectors and no rate lists, so the vehicle is uncontrolled.

## Design & Implementation
Every parameter is a keyword with a default except `planet`, `spacecraft`, `mission_time`, and `initial_time`, which are required keywords. Simulation settings disable plot generation and state normalisation while enabling results output to `_example_default_results_directory()`. Topography and wind are both switched off in the environment model. `solver_config` accepts `nothing`, letting the engine choose its own solver when the caller expresses no preference.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace re-exporting the SimulationModel configuration constructors. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `configuration` | SimulationConfiguration | n/a | — | Fully populated single-vehicle simulation configuration for examples and smoke runs. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:572-572`
- [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:612-612`
- [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callees` → `callers` · call · `src/precompile_workload.jl:22-22`

**Downstream**

- `callees` → [[analysis.example_support__example_default_results_directory|_example_default_results_directory]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:132-132`
- `callees` → [[core.simulation_configuration_environmentmodel|EnvironmentModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:151-151`
- `callees` → [[core.simulation_configuration_integrationtolerances|IntegrationTolerances]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:165-165`
- `callees` → [[core.simulation_configuration_missionconfiguration|MissionConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:143-143`
- `callees` → [[core.simulation_configuration_simulationsettings|SimulationSettings]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:136-136`
- `callees` → [[envana.env_gravity_models_inversesquaredj2gravitymodel|InverseSquaredJ2GravityModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:124-124`
- `callees` → [[environment.density_models_noatmospheremodel|NoAtmosphereModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:125-125`
- `callees` → [[environment.simple_ephemerides_spiceephemeridesmodel|SpiceEphemeridesModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:126-126`
- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:135-135`
- `callees` → [[vehicle.model_controlmodel|ControlModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:163-163`
- `callees` → [[vehicle.model_dynamicsmodel|DynamicsModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:160-160`
- `callees` → [[vehicle.model_guidancemodel|GuidanceModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:161-161`
- `callees` → [[vehicle.model_navigationmodel|NavigationModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:162-162`
- `callees` → [[vehicle.thermal_models_maxwellianheat|MaxwellianHeat]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:156-156`
<!-- vulcan:connections:end -->

## Limitations
Exactly one spacecraft is placed in the dynamics model, so multi-vehicle scenarios cannot be produced through this builder. The tolerance and step values are hard-coded rather than derived from the requested profile, so a strict verification run must build its configuration elsewhere. Because guidance, navigation, and control tuples are empty, any example needing closed-loop behaviour must construct its configuration directly.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/example_support.jl:119-175`, including the keyword defaults and the nested configuration constructors.
