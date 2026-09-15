---
id: environment.gravity_models_gravity_backbone_acceleration_ii
label: gravity_backbone_acceleration_ii
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: gravity_backbone_acceleration_ii
  lines:
  - 207
  - 207
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
  type: SVector{3,
  units: n/a
  description: Return value of `gravity_backbone_acceleration_ii`.
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

# gravity_backbone_acceleration_ii

## Purpose
Gravity-backbone acceleration hook (line 207) for `InverseSquaredGravityModel`: returns the inertial translational acceleration in m/s^2 that the backbone integrator applies in its gravity-only core. Overloads for `ConstantGravityModel` (identical) and `InverseSquaredJ2GravityModel` (planet-frame J2) follow in the file.

## Design & Implementation
Marked `@inline` with signature `(model, x::StateSample, env::EnvironmentSample, t::Float64)::SVector{3,Float64}`; simply returns `_inverse_squared_gravity_accel(x.pos_ii, env.planet)`. The J2 overload requires `env.planet_frame`, throwing `ArgumentError` otherwise, and computes `planet_frame.l_pi' * _inverse_squared_j2_gravity_accel(planet_frame.pos_pp, env.planet)` so the oblateness term is aligned with the true spin axis.

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
| out | `result` | SVector{3, | n/a | — | Return value of `gravity_backbone_acceleration_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:198-198`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:276-276`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`
- [[simulation.dynamics_rhs__gravity_backbone_core_acceleration|_gravity_backbone_core_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1502-1502`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:319-319`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/gravity/gravity_models.jl:218-218`
- `callees` → [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · feedback · `src/environment/gravity/gravity_models.jl:237-237`
- `callees` → [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/environment/gravity/gravity_models.jl:225-225`
- `callees` → [[dynamics.calc_force_torque|calcForceTorque]] · `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- `callees` → [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- `callees` → [[dynamics.perturbations_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/environment/gravity/gravity_models.jl:237-237`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/environment/gravity/gravity_models.jl:225-225`
- `callees` → [[dynamics.robot_arm_reaction_effector_calcforcetorque|calcForceTorque]] · `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- `callees` → [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callers` · call · `src/environment/gravity/gravity_models.jl:221-221`
- `callees` → [[environment.gravity_models__inverse_squared_gravity_accel|_inverse_squared_gravity_accel]] · `callers` · call · `src/environment/gravity/gravity_models.jl:213-213`
- `callees` → [[environment.gravity_models__inverse_squared_j2_gravity_accel|_inverse_squared_j2_gravity_accel]] · `callers` · call · `src/environment/gravity/gravity_models.jl:251-251`
- `callees` → [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callers` · call · `src/environment/gravity/gravity_models.jl:216-216`
- `callees` → [[environment.gravity_models_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/environment/gravity/gravity_models.jl:237-237`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/environment/gravity/gravity_models.jl:225-225`
<!-- vulcan:connections:end -->

## Limitations
Mass and torque are not part of the return, so this hook cannot express mass-dependent or rotational effects. `t` is unused, which is consistent with the position-only contract but means a time-varying `μ` is unsupported. For the J2 variant, `env.planet_frame` must have been sampled at the same stage time, and there is no check that `l_pi` corresponds to `t`.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 207.
