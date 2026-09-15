---
id: environment.gravity_models__inverse_squared_gravity_accel
label: _inverse_squared_gravity_accel
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: _inverse_squared_gravity_accel
  lines:
  - 38
  - 38
inputs:
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_inverse_squared_gravity_accel`.
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

# _inverse_squared_gravity_accel

## Purpose
Computes the point-mass gravitational acceleration vector at inertial position `pos_ii` for the given planet. It is the shared kernel behind `ConstantGravityModel`, `InverseSquaredGravityModel`, the aerobraking predictor's default gravity branch and the gravity-backbone acceleration hooks.

## Theory & Math
$$\mathbf a = -\frac{\mu}{r^2}\,\hat{\mathbf r} = -\frac{\mu}{|\mathbf r|^3}\,\mathbf r$$ with $\mu$ = `planet.μ` (m^3/s^2) and $\mathbf r$ = `pos_ii` (m).

## Design & Implementation
Marked `@inline`; takes `pos_ii::SVector{3,Float64}` in metres and any `planet` with a `μ` property. Computes `r = norm(pos_ii)`, `μ = Float64(planet.μ)`, and returns `-μ / r^2 * normalize(pos_ii)`, an `SVector{3,Float64}` in m/s^2. `normalize` divides by the norm a second time, so the expression is equivalent to `-μ * pos_ii / r^3`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_inverse_squared_gravity_accel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:163-163`
- [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:187-187`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:213-213`
- [[environment.gravity_models_wrench|wrench]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:199-199`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/gravity/gravity_models.jl:40-40`
<!-- vulcan:connections:end -->

## Limitations
At `r = 0` the result is `NaN` (0/0 inside `normalize`) with no guard; the gravity-gradient helper checks `r > 0` but this function does not. `norm` is evaluated twice (once explicitly, once inside `normalize`), a minor inefficiency on the hottest path in the simulator. No oblateness, no third-body terms, no relativistic corrections.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 38.
