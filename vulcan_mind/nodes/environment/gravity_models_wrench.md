---
id: environment.gravity_models_wrench
label: wrench
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: wrench
  lines:
  - 193
  - 193
inputs:
- id: model
  type: ConstantGravityModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: StateSample
  units: n/a
  required: true
  description: Positional argument `x`.
- id: env
  type: EnvironmentSample
  units: n/a
  required: true
  description: Positional argument `env`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `wrench`.
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

# wrench

## Purpose
Sampled-interface force/torque method (line 193) for `InverseSquaredGravityModel`, implementing the `EffectorSampling.wrench` hook. It returns the point-mass gravitational force in the inertial frame and the optional gravity-gradient torque in the body frame. Overloads for `ConstantGravityModel` (identical) and `InverseSquaredJ2GravityModel` (uses `env.planet_frame`) are defined in the same file.

## Design & Implementation
Marked `@inline`; signature `(model, x::StateSample, env::EnvironmentSample, t::Float64)`. Computes `gravity_ii = _inverse_squared_gravity_accel(x.pos_ii, env.planet)` and `force_ii = x.mass_kg * gravity_ii`, then `torque_body = _gravity_gradient_torque_body(model, x, env.planet)` which is zero unless `model.gravity_gradient` and both `x.q_ib` and `x.spacecraft` are present. The J2 overload throws `ArgumentError` when `env.planet_frame === nothing`, evaluates J2 on `planet_frame.pos_pp` and rotates back with `planet_frame.l_pi'`. `t` is unused.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | ConstantGravityModel | n/a | yes | Positional argument `model`. |
| in | `x` | StateSample | n/a | yes | Positional argument `x`. |
| in | `env` | EnvironmentSample | n/a | yes | Positional argument `env`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `wrench`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:154-154`
- [[core.effector_sampling_wrench_caching_bang|wrench_caching!]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:172-172`
- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1663-1663`
- [[dynamics.perturbations__lvlh_cascade_torque|_lvlh_cascade_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2263-2263`
- [[dynamics.perturbations__srp_total_acceleration_ii|_srp_total_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1365-1365`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2070-2070`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1925-1925`
- [[dynamics.perturbations_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1765-1765`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:259-259`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:225-225`
- [[grp.src_core_types|core/types/]] · `members_out` → `callers` · call · `src/core/types/effector_sampling.jl:154-154`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:274-274`

**Downstream**

- `callees` → [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callers` · call · `src/environment/gravity/gravity_models.jl:201-201`
- `callees` → [[environment.gravity_models__inverse_squared_gravity_accel|_inverse_squared_gravity_accel]] · `callers` · call · `src/environment/gravity/gravity_models.jl:199-199`
<!-- vulcan:connections:end -->

## Limitations
The return type annotation forces `Tuple{SVector{3,Float64}, SVector{3,Float64}}`, so any `planet.μ` stored as a non-`Float64` is converted inside the accel helper but a non-numeric value errors at that conversion. The gravity-gradient torque is computed from the central field only, even for the J2 model. Because the function is pure in `(model, x, env, t)`, mass loss during a step is not reflected until the next sample.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 193.
