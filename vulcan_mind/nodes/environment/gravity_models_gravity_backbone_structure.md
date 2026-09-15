---
id: environment.gravity_models_gravity_backbone_structure
label: gravity_backbone_structure
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: gravity_backbone_structure
  lines:
  - 205
  - 205
inputs:
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
  type: Symbol
  units: n/a
  description: Return value of `gravity_backbone_structure`. Returns `:position_only_static_gravity`.
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

# gravity_backbone_structure

## Purpose
Declares (line 205) that `InverseSquaredGravityModel` can participate in the gravity-backbone solver mode, whose translational core integrates only position-dependent static gravity. Equivalent one-line declarations exist for `ConstantGravityModel` and `InverseSquaredJ2GravityModel`.

## Design & Implementation
`@inline gravity_backbone_structure(::InverseSquaredGravityModel) = :position_only_static_gravity`, overriding the `:unsupported` fallback in `EffectorSampling`. The engine queries it once during solver setup and, on this answer, routes the model's acceleration through `gravity_backbone_acceleration_ii` instead of `wrench` for the backbone stage.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `gravity_backbone_structure`. Returns `:position_only_static_gravity`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_solver_partition|solver_partition]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:187-187`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:274-274`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:237-237`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`
- [[simulation.solver_policy__gravity_backbone_structure_validated|_gravity_backbone_structure_validated]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:150-150`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:307-307`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The declaration is unconditional: when `gravity_gradient = true` the torque part of the model is still not position-only in the rotational sense, and the backbone path returns acceleration only, so the gravity-gradient torque is dropped in that mode unless the engine separately calls `wrench` for torques. There is no mechanism to disable backbone participation per configuration.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 205.
