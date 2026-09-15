---
id: dynamics.perturbations_gravity_backbone_structure
label: gravity_backbone_structure
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: gravity_backbone_structure
  lines:
  - 1757
  - 1757
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
- dynamics
charts:
- dynamics
origin: agent
---

# gravity_backbone_structure

## Purpose
Tells the split solver that spherical-harmonics gravity is a position-only static field, qualifying it for the implicit backbone that handles the stiff central attraction.

## Design & Implementation
Returns the symbol `:position_only_static_gravity`. Declared `@inline`. Position-only means the acceleration depends on position but not velocity or time, which the implicit partition's Jacobian assumptions require.

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
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.solver_policy__gravity_backbone_structure_validated|_gravity_backbone_structure_validated]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:150-150`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:307-307`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Strictly, the harmonics field is time-dependent through the planet's rotation; the classification treats that dependence as slow enough to ignore within one step.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1757.
