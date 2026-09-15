---
id: environment.density_models_polynomialfitatmospheremodel
label: PolynomialFitAtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: PolynomialFitAtmosphereModel
  lines:
  - 199
  - 199
inputs:
- id: polyfit_coeffs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `polyfit_coeffs`.
- id: valid_min_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Field `valid_min_altitude_m`.
- id: valid_max_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Field `valid_max_altitude_m`.
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
  type: PolynomialFitAtmosphereModel
  units: n/a
  description: Constructed `PolynomialFitAtmosphereModel`.
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

# PolynomialFitAtmosphereModel

## Purpose
An analytic model fitting log-density as a polynomial in altitude in kilometres, used for the smooth above-interface region where GRAM sampling would be wasted.

## Design & Implementation
Immutable with a coefficient vector in Horner order and validity bounds that are enforced as a clamp on the evaluation altitude. Log-density is further clamped to the finite `exp` range of `Float64`. A planet constructor reads the coefficients and optional bounds from the planet object, with defaults of 50 km to 2,000 km. Temperature comes from the planet's `T_ref`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `polyfit_coeffs` | Vector{Float64} | n/a | yes | Field `polyfit_coeffs`. |
| in | `valid_min_altitude_m` | Float64 | n/a | yes | Field `valid_min_altitude_m`. |
| in | `valid_max_altitude_m` | Float64 | n/a | yes | Field `valid_max_altitude_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PolynomialFitAtmosphereModel | n/a | — | Constructed `PolynomialFitAtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__planet_polyfit_valid_max_altitude_m|_planet_polyfit_valid_max_altitude_m]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:232-232`
- [[environment.density_models_density_polyfit|density_polyfit]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1082-1082`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:210-210`
<!-- vulcan:connections:end -->

## Limitations
Clamping at the bounds makes density constant beyond them rather than decaying, which is physically wrong but numerically safe; the polynomial has no latitude, longitude or time dependence.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 199.
