---
id: environment.density_models__gram_point_density
label: _gram_point_density
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _gram_point_density
  lines:
  - 1036
  - 1036
inputs:
- id: model
  type: Any
  units: n/a
  required: true
  description: Positional argument `model`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
- id: lat
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat`.
- id: lon
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon`.
- id: el_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `el_time`.
- id: wind
  type: Bool
  units: n/a
  required: true
  description: Positional argument `wind`.
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_gram_point_density`.
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

# _gram_point_density

## Purpose
Evaluates the native GRAM point model at one location, the fallback the surrogate uses below its point-fallback altitude.

## Design & Implementation
The core package provides a catch-all method that throws `MethodError` for unknown model types and a forwarding method for `GRAMAtmosphereModelSurrogate` that delegates to its base model. The extension adds the real method for `GRAMAtmosphereModel`, clamping altitude at −30 m and calling `GRAMSuite.point_density_state` under the model's lock; the telemetry fallback base type adds a vacuum-returning method.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | Any | n/a | yes | Positional argument `model`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `lat` | Float64 | n/a | yes | Positional argument `lat`. |
| in | `lon` | Float64 | n/a | yes | Positional argument `lon`. |
| in | `el_time` | Float64 | n/a | yes | Positional argument `el_time`. |
| in | `wind` | Bool | n/a | yes | Positional argument `wind`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_gram_point_density`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__gramofflinesurrogatefallbackbase|_GRAMOfflineSurrogateFallbackBase]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:277-277`
- [[ext.gram_core_density_state|_gram_core_density_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:255-255`
- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:255-255`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_density_state|_nrlmsise_density_state]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1078-1078`
- `callees` → [[environment.density_models__nrlmsise_eval_datetime|_nrlmsise_eval_datetime]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1077-1077`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/environment/atmosphere/density_models.jl:1058-1058`
<!-- vulcan:connections:end -->

## Limitations
The explicit `MethodError` throw in the catch-all hides which types are supported from the method table; a reader must know to look in the extension.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 1036.
