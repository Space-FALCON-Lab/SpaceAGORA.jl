---
id: environment.density_models__piecewise_layer_index
label: _piecewise_layer_index
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _piecewise_layer_index
  lines:
  - 478
  - 478
inputs:
- id: model
  type: PiecewiseExponentialAtmosphereModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
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
  type: Int
  units: n/a
  description: Return value of `_piecewise_layer_index`.
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

# _piecewise_layer_index

## Purpose
Selects which layer of a piecewise exponential model applies at an altitude, extrapolating the end layers outside the breakpoints.

## Design & Implementation
Uses `searchsortedlast` on the breakpoints and clamps the index into one through the layer count. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | PiecewiseExponentialAtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_piecewise_layer_index`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:962-962`
- [[environment.density_models_timetabulatedatmospheremodel|TimeTabulatedAtmosphereModel]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:824-824`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Binary search per evaluation; for the handful of layers typically used a linear scan would be equally fast, but the search keeps large tables efficient.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 478.
