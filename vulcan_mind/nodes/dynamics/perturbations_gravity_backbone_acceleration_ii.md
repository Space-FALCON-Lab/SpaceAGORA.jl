---
id: dynamics.perturbations_gravity_backbone_acceleration_ii
label: gravity_backbone_acceleration_ii
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: gravity_backbone_acceleration_ii
  lines:
  - 1759
  - 1759
inputs:
- id: model
  type: GravitationalHarmonicsModel
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
- dynamics
charts:
- dynamics
origin: agent
---

# gravity_backbone_acceleration_ii

## Purpose
Exposes the spherical-harmonics acceleration through the gravity-backbone interface that the split and implicit solver partitions use for the position-only static field.

## Design & Implementation
Calls `wrench(model, x, env, t)`, discards the torque, and divides the force by `x.mass_kg`. Declared `@inline`. The backbone interface deals in accelerations because the implicit partition integrates the state without mass.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
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
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__gravity_backbone_core_acceleration|_gravity_backbone_core_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1502-1502`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:319-319`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1765-1765`
- `callees` → [[dynamics.perturbations_get_magnetic_field_dipole|get_magnetic_field_dipole]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1770-1770`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1765-1765`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1765-1765`
<!-- vulcan:connections:end -->

## Limitations
Ignores the torque component, which is zero for this model in any case; and pays the full wrench evaluation even though only force is needed.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1759.
