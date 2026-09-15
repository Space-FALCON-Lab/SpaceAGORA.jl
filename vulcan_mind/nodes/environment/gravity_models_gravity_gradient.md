---
id: environment.gravity_models_gravity_gradient
label: gravity_gradient
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: gravity_gradient
  lines:
  - 303
  - 303
inputs:
- id: J
  type: SMatrix{3,3,Float64}
  units: n/a
  required: true
  description: Positional argument `J`.
- id: rVec
  type: SVector{3,Float64}
  units: n/a
  required: true
  description: Positional argument `rVec`.
- id: mu
  type: Float64
  units: n/a
  required: true
  description: Positional argument `μ`.
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
  type: SVector
  units: n/a
  description: Return value of `gravity_gradient`. Returns `SVector{3, Float64}(0.0,
    0.0, 0.0)` or `3*μ/r^3 * cross(r_hat, J * r_hat)`.
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

# gravity_gradient

## Purpose
Computes the classical gravity-gradient torque acting on a rigid body in a central gravitational field, expressed in the body frame, given the body-frame inertia tensor and position vector. It is the physics kernel called by both `_gravity_gradient_torque_body` methods.

## Theory & Math
$$\boldsymbol\tau_{gg} = \frac{3\mu}{r^3}\,\hat{\mathbf r} \times \left(\mathbf J\,\hat{\mathbf r}\right)$$ where $\mathbf J$ is the body-frame inertia tensor (kg m^2), $\hat{\mathbf r}$ the unit position vector in the body frame, $r$ its magnitude (m) and $\mu$ the gravitational parameter (m^3/s^2).

## Design & Implementation
Signature `gravity_gradient(J::SMatrix{3,3,Float64}, rVec::SVector{3,Float64}, μ::Float64)`. Computes `r = norm(rVec)` and returns the zero vector if `r` is non-finite or non-positive. Otherwise forms `r_hat = rVec / r` and returns `3μ/r^3 * cross(r_hat, J * r_hat)` as an `SVector{3,Float64}` in N m. Being allocation-free with static arrays, it is safe on the ODE hot path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `J` | SMatrix{3,3,Float64} | n/a | yes | Positional argument `J`. |
| in | `rVec` | SVector{3,Float64} | n/a | yes | Positional argument `rVec`. |
| in | `mu` | Float64 | n/a | yes | Positional argument `μ`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `gravity_gradient`. Returns `SVector{3, Float64}(0.0, 0.0, 0.0)` or `3*μ/r^3 * cross(r_hat, J * r_hat)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:89-89`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:289-289`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Valid only for a point-mass central field; oblateness contributions to the gradient torque are omitted even when the force model includes J2. The formula assumes the spacecraft dimensions are small compared with `r` (first-order expansion of the field). No return-type annotation is given, so a non-`SMatrix` inertia argument would fail to dispatch rather than convert. The zero-return for degenerate `r` silently suppresses errors.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 303.
