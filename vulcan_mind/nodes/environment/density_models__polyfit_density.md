---
id: environment.density_models__polyfit_density
label: _polyfit_density
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _polyfit_density
  lines:
  - 498
  - 498
inputs:
- id: model
  type: PolynomialFitAtmosphereModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: h_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h_m`.
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
  type: Float64
  units: n/a
  description: Return value of `_polyfit_density`.
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

# _polyfit_density

## Purpose
Returns the polynomial-fit model's density in kilograms per cubic metre at a given altitude.

## Design & Implementation
Exponentiates the result of `_polyfit_log_density`, which has already clamped both the altitude to the validity band and the log-density to the finite range of `exp` for `Float64`. Declared `@inline`. Used by both the scalar `getDensity` and the batch method for `PolynomialFitAtmosphereModel`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | PolynomialFitAtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `h_m` | Float64 | n/a | yes | Positional argument `h_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_polyfit_density`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:986-986`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__polyfit_log_density|_polyfit_log_density]] · `callers` · call · `src/environment/atmosphere/density_models.jl:499-499`
<!-- vulcan:connections:end -->

## Limitations
Inherits the clamped-band behaviour, and because the log-density floor is `log(nextfloat(0.0))` the smallest representable density is a subnormal rather than exactly zero, so drag never vanishes entirely under this model.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 498.
