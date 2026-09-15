---
id: environment.gravity_models__inverse_squared_j2_gravity_accel
label: _inverse_squared_j2_gravity_accel
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: _inverse_squared_j2_gravity_accel
  lines:
  - 44
  - 44
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
  description: Return value of `_inverse_squared_j2_gravity_accel`.
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

# _inverse_squared_j2_gravity_accel

## Purpose
Computes gravitational acceleration including the J2 zonal oblateness term, for use by `InverseSquaredJ2GravityModel` and by the aerobraking predictor when `gm_code == 2`. The input position must be expressed in a frame whose z axis is the planet's spin axis (the planet-fixed frame), which is why the `wrench` method feeds it `planet_frame.pos_pp`.

## Theory & Math
$$\mathbf a = -\frac{\mu}{r^2}\hat{\mathbf r} + \frac{3}{2} J_2 \frac{\mu R_e^2}{r^4} \begin{pmatrix} \frac{x}{r}\left(5\frac{z^2}{r^2} - 1\right) \\ \frac{y}{r}\left(5\frac{z^2}{r^2} - 1\right) \\ \frac{z}{r}\left(5\frac{z^2}{r^2} - 3\right) \end{pmatrix}$$ where $R_e$ = `planet.Rp_e`, $J_2$ = `planet.J2`, and $(x, y, z)$ are body-fixed coordinates.

## Design & Implementation
Marked `@inline`. Reads `μ`, `J2` and the equatorial radius `Rp_e` from `planet` (the comment notes that J2 is normalised to the equatorial radius to match `j2_secular_rates` and the Pines harmonics path). With `x, y, z = pos_ii`, `r^2` and `r` it forms the spherical term `-μ/r^2 * r_hat` and the J2 vector `(x/r (5z^2/r^2 - 1), y/r (5z^2/r^2 - 1), z/r (5z^2/r^2 - 3))` scaled by `3/2 J2 μ Rp^2 / r^4`, returning their sum as `SVector{3,Float64}` in m/s^2.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_inverse_squared_j2_gravity_accel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:160-160`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:267-267`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:251-251`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/gravity/gravity_models.jl:46-46`
<!-- vulcan:connections:end -->

## Limitations
The signature names the argument `pos_ii` and the legacy `calcForceTorque` path passes an inertial position, so in that path the J2 term is oriented along the inertial z axis rather than the true spin axis; only the `wrench` and backbone paths rotate correctly. There is no guard for `r = 0`. Only J2 is modelled; J3, J4 and tesseral terms require the Pines harmonics path. `Rp_e` being the equatorial radius must be consistent with how `planet.J2` was derived.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 44.
