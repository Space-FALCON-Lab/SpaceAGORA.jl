---
id: environment.density_models__density_scalar_for_batch
label: _density_scalar_for_batch
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _density_scalar_for_batch
  lines:
  - 864
  - 864
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
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_density_scalar_for_batch`.
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

# _density_scalar_for_batch

## Purpose
Routes each element of a generic batch density call to the correct scalar `getDensity` signature for the model type.

## Design & Implementation
Three `@inline` methods: the generic one forwards all seven arguments including `p`; the `ExponentialAtmosphereModel` method drops `p` because that model has no parameter-dependent form; the `NRLMSISE00AtmosphereModel` method keeps `p` because the epoch is required.

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
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_density_scalar_for_batch`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1011-1011`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/environment/atmosphere/density_models.jl:873-873`
<!-- vulcan:connections:end -->

## Limitations
The generic batch path is a serial loop over scalar calls, so models with expensive per-call setup gain nothing from batching unless they specialise `getDensityBatch!` themselves.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 864.
