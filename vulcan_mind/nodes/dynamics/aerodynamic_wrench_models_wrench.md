---
id: dynamics.aerodynamic_wrench_models_wrench
label: wrench
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: wrench
  lines:
  - 512
  - 512
inputs:
- id: model
  type: AerodynamicCoefficientConstant
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
- dynamics
charts:
- dynamics
origin: agent
---

# wrench

## Purpose
Stateless effector interface returning inertial force and root-body-frame torque for the three aerodynamic models from a `StateSample` and `EnvironmentSample`.

## Design & Implementation
Three `@inline` methods with signature `(model, x::StateSample, env::EnvironmentSample, t::Float64)`. `AerodynamicCoefficientConstant` and `AerodynamicCoefficientNoBallisticFlight` call `_aero_pure_wrench(:constant, x, env)`; `AerodynamicCoefficientfM` calls `_aero_pure_wrench(:fm, x, env, nothing, model.fixed_attitude_incidence)`. Each discards the drag/lift/cross components and returns `(force, torque)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerodynamicCoefficientConstant | n/a | yes | Positional argument `model`. |
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

- `callees` → [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:518-518`
<!-- vulcan:connections:end -->

## Limitations
Without `ODEParams` these methods cannot perform per-link atmosphere sampling, so `per_link_atmosphere=true` has no effect on this path. `t` is unused. The `NoBallisticFlight` variant is a duplicate of the constant model rather than Newtonian flow.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 512.
