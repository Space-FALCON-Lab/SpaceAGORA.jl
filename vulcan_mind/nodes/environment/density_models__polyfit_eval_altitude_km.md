---
id: environment.density_models__polyfit_eval_altitude_km
label: _polyfit_eval_altitude_km
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _polyfit_eval_altitude_km
  lines:
  - 483
  - 483
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
  description: Return value of `_polyfit_eval_altitude_km`.
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

# _polyfit_eval_altitude_km

## Purpose
Clamps a query altitude into the polynomial model's validity band and converts it to kilometres, the unit the polynomial coefficients were fitted in.

## Design & Implementation
Returns `clamp(h_m, valid_min_altitude_m, valid_max_altitude_m) * 1e-3`. Declared `@inline`. The clamp is deliberate: an unbounded polynomial in altitude diverges rapidly outside its fit domain, so evaluation is pinned to the band edges rather than extrapolated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | PolynomialFitAtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `h_m` | Float64 | n/a | yes | Positional argument `h_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_polyfit_eval_altitude_km`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__polyfit_log_density|_polyfit_log_density]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:490-490`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Pinning makes density constant outside the band instead of decaying, which is physically wrong above the upper bound but numerically safe; a caller that needs correct high-altitude behaviour must choose a different model.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 483.
