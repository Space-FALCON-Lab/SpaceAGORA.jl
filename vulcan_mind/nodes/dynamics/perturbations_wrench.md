---
id: dynamics.perturbations_wrench
label: wrench
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: wrench
  lines:
  - 993
  - 993
inputs:
- id: model
  type: NBodyGravityModel
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
The generic wrench-interface implementation for the N-body effector, returning force from the precomputed third-body positions in the environment sample and zero torque.

## Design & Implementation
Requires `env.third_bodies` to be present, raising `ArgumentError` otherwise, computes the acceleration with `_nbody_acceleration_ii` and multiplies by `x.mass_kg`. Returns the force and a zero torque vector. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | NBodyGravityModel | n/a | yes | Positional argument `model`. |
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

- `callees` → [[dynamics.perturbations__nbody_acceleration_ii|_nbody_acceleration_ii]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1001-1001`
<!-- vulcan:connections:end -->

## Limitations
Raises when the environment sample lacks third bodies, which happens only if the sampler ignored this effector's `environment_requirements`; there is no fallback to a live SPICE lookup here.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 993.
