---
id: environment.density_models_gramatmospheremodelsurrogate
label: GRAMAtmosphereModelSurrogate
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: GRAMAtmosphereModelSurrogate
  lines:
  - 179
  - 179
inputs:
- id: base_model
  type: M
  units: n/a
  required: true
  description: Field `base_model`.
- id: surrogate_file
  type: String
  units: n/a
  required: true
  description: Field `surrogate_file`.
- id: point_fallback_below_m
  type: Union{Nothing, Float64}
  units: n/a
  required: true
  description: Field `point_fallback_below_m`.
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
  type: GRAMAtmosphereModelSurrogate
  units: n/a
  description: Constructed `GRAMAtmosphereModelSurrogate`.
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

# GRAMAtmosphereModelSurrogate

## Purpose
A GRAM surrogate wrapper that answers density from a precomputed offline grid file, falling back to a point model below a configurable altitude.

## Design & Implementation
An immutable struct parameterised on the base model type, holding `base_model`, the `surrogate_file` path and an optional `point_fallback_below_m` altitude. The base is normally a `GRAMAtmosphereModel` but may be any type implementing `_gram_point_density`, which is how the library-less telemetry fallback works. Property access forwards to the base model. The extension supplies the constructor and `getDensity`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `base_model` | M | n/a | yes | Field `base_model`. |
| in | `surrogate_file` | String | n/a | yes | Field `surrogate_file`. |
| in | `point_fallback_below_m` | Union{Nothing, Float64} | n/a | yes | Field `point_fallback_below_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GRAMAtmosphereModelSurrogate | n/a | — | Constructed `GRAMAtmosphereModelSurrogate`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__try_libraryless_gram_surrogate|_try_libraryless_gram_surrogate]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:313-313`
- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:126-126`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The surrogate grid is frozen at one epoch and cannot represent the diurnal or dust-storm variability of the native model; the point fallback altitude is the only knob for reintroducing native fidelity near periapsis.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 179.
