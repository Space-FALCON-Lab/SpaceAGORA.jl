---
id: core.effector_sampling_gravity_backbone_structure
label: gravity_backbone_structure
kind: function
source:
  file: src/core/types/effector_sampling.jl
  symbol: gravity_backbone_structure
  lines:
  - 195
  - 195
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
  description: Return value of `gravity_backbone_structure`. Returns `:unsupported`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# gravity_backbone_structure

## Purpose
Declaration hook for the gravity-backbone solver mode. An effector returns `:position_only_static_gravity` to advertise that its acceleration depends only on position and static planet data, allowing it to be folded into the fast gravity-only translational core; otherwise it stays on the general path.

## Design & Implementation
The fallback `@inline gravity_backbone_structure(::Any) = :unsupported` opts every model out. A model that opts in must also implement `gravity_backbone_acceleration_ii(model, x::StateSample, env::EnvironmentSample, t)` returning inertial acceleration in m/s^2. The engine queries this once when selecting the backbone participants.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `gravity_backbone_structure`. Returns `:unsupported`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_solver_partition|solver_partition]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:187-187`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:274-274`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · feedback · `src/environment/gravity/gravity_models.jl:237-237`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`
- [[simulation.solver_policy__gravity_backbone_structure_validated|_gravity_backbone_structure_validated]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:150-150`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:307-307`

**Downstream**

- `callees` → [[core.effector_sampling_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `callers` · call · `src/core/types/effector_sampling.jl:209-209`
- `callees` → [[dynamics.perturbations_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/core/types/effector_sampling.jl:198-198`
- `callees` → [[dynamics.perturbations_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `callers` · call · `src/core/types/effector_sampling.jl:209-209`
- `callees` → [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callers` · call · `src/core/types/effector_sampling.jl:198-198`
<!-- vulcan:connections:end -->

## Limitations
Nothing verifies that an opted-in effector really is position-only; a model that reads `x.vel_ii` inside `gravity_backbone_acceleration_ii` would silently violate the symplectic assumptions of the backbone integrator. The missing-method failure for a model that declares support but lacks the acceleration hook occurs at the first ODE evaluation, not at configuration time. Only two symbols are meaningful and neither is validated.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 195.
