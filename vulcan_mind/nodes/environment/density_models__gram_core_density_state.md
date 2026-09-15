---
id: environment.density_models__gram_core_density_state
label: _gram_core_density_state
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _gram_core_density_state
  lines:
  - 461
  - 461
inputs:
- id: _core
  type: Any
  units: n/a
  required: true
  description: Positional argument `_core`.
- id: _h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `_h`.
- id: _lat
  type: Float64
  units: n/a
  required: true
  description: Positional argument `_lat`.
- id: _lon
  type: Float64
  units: n/a
  required: true
  description: Positional argument `_lon`.
- id: _el_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `_el_time`.
- id: _wind
  type: Bool
  units: n/a
  required: true
  description: Positional argument `_wind`.
- id: _lock_obj
  type: Any
  units: n/a
  required: true
  description: Positional argument `_lock_obj`.
- id: _vacuum_temperature
  type: Float64
  units: n/a
  required: true
  description: Positional argument `_vacuum_temperature`.
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
  description: Return value of `_gram_core_density_state`.
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

# _gram_core_density_state

## Purpose
The core-package stub for evaluating a native GRAM core, overridden by the extension's real method for `GRAMSuite.GRAMAtmosphereModel` cores.

## Design & Implementation
Accepts the core, altitude, latitude, longitude, elapsed time, wind flag, lock object and vacuum temperature, and raises the not-loaded error. The isolated-pool batch path calls this signature so that pool workers can evaluate a core directly with a chosen lock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `_core` | Any | n/a | yes | Positional argument `_core`. |
| in | `_h` | Float64 | n/a | yes | Positional argument `_h`. |
| in | `_lat` | Float64 | n/a | yes | Positional argument `_lat`. |
| in | `_lon` | Float64 | n/a | yes | Positional argument `_lon`. |
| in | `_el_time` | Float64 | n/a | yes | Positional argument `_el_time`. |
| in | `_wind` | Bool | n/a | yes | Positional argument `_wind`. |
| in | `_lock_obj` | Any | n/a | yes | Positional argument `_lock_obj`. |
| in | `_vacuum_temperature` | Float64 | n/a | yes | Positional argument `_vacuum_temperature`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_gram_core_density_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:234-234`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`
- [[simulation.model_selection__gram_isolated_pool_density_state|_gram_isolated_pool_density_state]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:103-103`

**Downstream**

- `callees` → [[environment.density_models__gram_not_loaded_error|_gram_not_loaded_error]] · `callers` · call · `src/environment/atmosphere/density_models.jl:471-471`
<!-- vulcan:connections:end -->

## Limitations
Every argument is prefixed with an underscore to silence unused-variable warnings, which makes the stub's intended signature harder to read than the extension's method.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 461.
