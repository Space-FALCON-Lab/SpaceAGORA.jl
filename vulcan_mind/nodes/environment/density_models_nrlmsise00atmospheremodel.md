---
id: environment.density_models_nrlmsise00atmospheremodel
label: NRLMSISE00AtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: NRLMSISE00AtmosphereModel
  lines:
  - 309
  - 309
inputs:
- id: f107a
  type: Float64
  units: n/a
  required: true
  description: Field `f107a`.
- id: f107
  type: Float64
  units: n/a
  required: true
  description: Field `f107`.
- id: ap
  type: A
  units: n/a
  required: true
  description: Field `ap`.
- id: index_provider
  type: P
  units: n/a
  required: true
  description: Field `index_provider`.
- id: include_anomalous_oxygen
  type: Bool
  units: n/a
  required: true
  description: Field `include_anomalous_oxygen`.
- id: valid_min_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Field `valid_min_altitude_m`.
- id: valid_max_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Field `valid_max_altitude_m`.
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
  type: NRLMSISE00AtmosphereModel
  units: n/a
  description: Constructed `NRLMSISE00AtmosphereModel`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# NRLMSISE00AtmosphereModel

## Purpose
Earth's NRLMSISE-00 empirical atmosphere with either fixed geophysical indices or a callable provider, wrapping `SatelliteToolboxAtmosphericModels`.

## Design & Implementation
Immutable, parameterised on the Ap type (scalar or seven-vector) and provider type, with `f107a`, `f107`, `ap`, `index_provider`, the anomalous-oxygen flag and validity bounds. The keyword constructor rejects combining `use_space_indices` with a custom provider, requires the download flag only with space indices, validates fluxes and Ap, and installs `NRLMSISE00SpaceIndicesProvider` when requested. The six-argument `getDensity` deliberately throws because elapsed time without an epoch would anchor to J2000; the seven-argument form maps elapsed time through the run's `initial_time`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f107a` | Float64 | n/a | yes | Field `f107a`. |
| in | `f107` | Float64 | n/a | yes | Field `f107`. |
| in | `ap` | A | n/a | yes | Field `ap`. |
| in | `index_provider` | P | n/a | yes | Field `index_provider`. |
| in | `include_anomalous_oxygen` | Bool | n/a | yes | Field `include_anomalous_oxygen`. |
| in | `valid_min_altitude_m` | Float64 | n/a | yes | Field `valid_min_altitude_m`. |
| in | `valid_max_altitude_m` | Float64 | n/a | yes | Field `valid_max_altitude_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NRLMSISE00AtmosphereModel | n/a | — | Constructed `NRLMSISE00AtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:189-189`
- [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:342-342`
- [[environment.density_models__planet_polyfit_valid_max_altitude_m|_planet_polyfit_valid_max_altitude_m]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:255-255`
- [[environment.density_models_nrlmsise00spaceindicesprovider|NRLMSISE00SpaceIndicesProvider]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:272-272`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:218-218`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Earth only, and winds are always zero; validity bounds are documentation, not a clamp, so evaluation well above 1,000 km returns whatever the model extrapolates.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 309.
