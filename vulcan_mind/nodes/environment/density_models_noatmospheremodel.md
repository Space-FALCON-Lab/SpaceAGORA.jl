---
id: environment.density_models_noatmospheremodel
label: NoAtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: NoAtmosphereModel
  lines:
  - 12
  - 12
inputs:
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
  type: NoAtmosphereModel
  units: n/a
  description: Constructed `NoAtmosphereModel`.
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

# NoAtmosphereModel

## Purpose
The vacuum density model: zero density everywhere, used by the no-GRAM quickstart and by any run that wants drag switched off without removing the aerodynamic effector.

## Design & Implementation
An empty immutable subtype of `AbstractDensityModel`. Its `getDensity` returns zero density, the planet's `T_ref` and zero wind, and its batch method fills the output vectors likewise. It is the one model for which `density_vanishes_above_entry_interface` returns true, so the split solver may skip aerodynamic evaluation on coast arcs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NoAtmosphereModel | n/a | — | Constructed `NoAtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:177-177`
- [[core.no_gram_presets_make_no_gram_density_model|make_no_gram_density_model]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:53-53`
- [[core.no_gram_presets_make_no_gram_planet|make_no_gram_planet]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:41-41`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:125-125`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:187-187`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because it reports `T_ref` rather than a physical temperature, thermal models fed from it see a constant reference value that has no relationship to any real exosphere.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 12.
