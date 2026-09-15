---
id: environment.density_models__polyfit_log_density
label: _polyfit_log_density
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _polyfit_log_density
  lines:
  - 487
  - 487
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
  description: Return value of `_polyfit_log_density`.
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

# _polyfit_log_density

## Purpose
Evaluates the log-density polynomial by Horner's rule and bounds the result to the finite exponential range.

## Theory & Math
$$
\ln\rho = \sum_{k=0}^{n} c_k\, h_{\text{km}}^{\,n-k}\quad\text{evaluated as}\quad ((c_0 h + c_1) h + c_2) h + \cdots
$$

## Design & Implementation
Returns zero for an empty coefficient vector. Otherwise it starts from the leading coefficient and folds the rest with `muladd` against the clamped altitude in kilometres, then clamps into `[log(nextfloat(0)), log(floatmax)]`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | PolynomialFitAtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `h_m` | Float64 | n/a | yes | Positional argument `h_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_polyfit_log_density`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__polyfit_density|_polyfit_density]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:499-499`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__polyfit_eval_altitude_km|_polyfit_eval_altitude_km]] · `callers` · call · `src/environment/atmosphere/density_models.jl:490-490`
<!-- vulcan:connections:end -->

## Limitations
An empty coefficient vector yields `exp(0) = 1 kg/m³` everywhere, which is a surprising sentinel rather than an error.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 487.
