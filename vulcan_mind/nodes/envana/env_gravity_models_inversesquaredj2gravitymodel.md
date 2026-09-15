---
id: envana.env_gravity_models_inversesquaredj2gravitymodel
label: InverseSquaredJ2GravityModel
kind: struct
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: InverseSquaredJ2GravityModel
  lines:
  - 18
  - 23
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Environment namespace providing the J2 acceleration and gravity-gradient
    torque kernels.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: wrench
  type: Tuple{SVector{3,Float64},SVector{3,Float64}}
  units: N,N*m
  description: Inertial gravity force and body-frame gravity-gradient torque for one
    vehicle.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- envana
origin: agent
---
# InverseSquaredJ2GravityModel

## Purpose
`InverseSquaredJ2GravityModel` is the oblate-planet gravity effector type. Declaring it selects point-mass gravity augmented with the second zonal harmonic, and its single field decides whether the accompanying gravity-gradient torque is produced in addition to the force.

## Theory & Math
The acceleration combines the point-mass and second-zonal terms. Writing `r` for the inertial position in metres, `mu` for the gravitational parameter in m^3/s^2, `R_e` for the equatorial radius in metres, and `J2` for the dimensionless oblateness coefficient, the acceleration is `a = -mu r / |r|^3 + (1.5 J2 mu R_e^2 / |r|^5) * v`, where `v` collects the `(5 (z/|r|)^2 - 1)` factors on the equatorial components and `(5 (z/|r|)^2 - 3)` on the axial one. Force follows as `F = m a` in newtons with mass `m` in kilograms. The gradient torque is `T = 3 mu / |r|^3 * (u_hat x (J u_hat))`, with inertia tensor `J` in kg*m^2 and `u_hat` the unit nadir direction expressed in the body frame; its units are N*m.

## Model & Assumptions
The struct carries no physical constants of its own: `mu`, `R_e`, and `J2` are read from `param.args.environment_model.planet` at evaluation time, and the commented-out Earth values at lines 19 through 21 record the magnitudes involved, namely `3.986004418e14` m^3/s^2, `1.08263e-3`, and `6378137.0` m. The gradient torque assumes a rigid body with a constant inertia tensor and the small-body approximation in which the vehicle dimension is negligible compared with `|r|`. It is emitted only when `gravity_gradient` is true and the mission configuration has `orientation_sim` enabled.

## Design & Implementation
The type is declared with `@kwdef` and subtypes `AbstractForceTorqueModel`, so it is constructed with no arguments in the common case and dispatches through the standard effector interface. Its `calcForceTorque` method at line 248 unpacks position from state entries 1 through 3 and mass from entry 7, then delegates to `_inverse_squared_j2_gravity_accel` for the acceleration and `_gravity_gradient_torque_body` for the torque, both returning `SVector{3,Float64}` so the whole evaluation is allocation-free inside the ODE right-hand side. A companion `environment_requirements` method at line 257 declares `planet_frame=true`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Environment namespace providing the J2 acceleration and gravity-gradient torque kernels. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `wrench` | Tuple{SVector{3,Float64},SVector{3,Float64}} | N,N*m | — | Inertial gravity force and body-frame gravity-gradient torque for one vehicle. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__base_gravity_effector|_base_gravity_effector]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:24-24`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:124-124`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Gravity is truncated at `J2`; higher zonal and all tesseral terms are absent, so long-arc ground-track prediction drifts. The torque model uses the point-nadir approximation and ignores the gradient of `J2` itself. Only the vehicle at index `i` is evaluated, so multi-body mutual gravitation is out of scope for this method.

## Provenance
Read directly from `src/environment/gravity/gravity_models.jl:18-23`, together with its `calcForceTorque` dispatch at line 248 and the shared `gravity_gradient` kernel at line 303 in the same file.
