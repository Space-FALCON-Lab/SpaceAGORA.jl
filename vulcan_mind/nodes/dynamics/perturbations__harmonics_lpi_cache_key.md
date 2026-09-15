---
id: dynamics.perturbations__harmonics_lpi_cache_key
label: _harmonics_lpi_cache_key
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _harmonics_lpi_cache_key
  lines:
  - 560
  - 560
inputs:
- id: model
  type: GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: param
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `param`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
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
  type: Any
  units: n/a
  description: Return value of `_harmonics_lpi_cache_key`. Returns `(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _harmonics_lpi_cache_key

## Purpose
Forms the tuple that identifies one planet-frame rotation request, so the one-entry `harmonics_lpi` cache on shared buffers can tell whether the stored matrix is the one being asked for.

## Design & Implementation
Returns `(model.planet.name, ephemerides_cache_key(ephemerides_model), et)`. Including the ephemerides model key means a run that switches between SPICE and analytic frame models cannot be served a matrix from the other; including the planet name guards multi-planet configurations. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_harmonics_lpi_cache_key`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__harmonics_lpi_at_bang|_harmonics_lpi_at!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:574-574`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[environment.simple_ephemerides_ephemerides_cache_key|ephemerides_cache_key]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:563-563`
<!-- vulcan:connections:end -->

## Limitations
The comparison on `et` is exact floating-point equality, so RK stage times that differ from the cached one by a single ulp miss the cache and pay a fresh `planet_frame_lpi` under the lock.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 560.
