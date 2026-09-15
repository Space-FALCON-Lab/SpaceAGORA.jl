---
id: core.effector_sampling_gravity_backbone_kick_structure
label: gravity_backbone_kick_structure
kind: function
source:
  file: src/core/types/effector_sampling.jl
  symbol: gravity_backbone_kick_structure
  lines:
  - 218
  - 218
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
  description: Return value of `gravity_backbone_kick_structure`. Returns `:unsupported`.
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

# gravity_backbone_kick_structure

## Purpose
Declaration hook for the `gravity_backbone_split` solver mode. It lets a translational perturbation effector announce that it should be applied as an explicit velocity kick around the gravity core rather than being integrated inside it, which is the operator-splitting treatment for weak, smooth perturbations such as SRP or third-body gravity.

## Design & Implementation
Default `@inline gravity_backbone_kick_structure(::Any) = :unsupported`. Effectors overload it to return `:velocity_kick_explicit` and must then implement `gravity_backbone_kick_acceleration_ii(model, x, env, t)` returning inertial acceleration in m/s^2, which the engine multiplies by the substep to form the kick. The hook is evaluated once during solver setup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `gravity_backbone_kick_structure`. Returns `:unsupported`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:209-209`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`
- [[simulation.solver_policy__gravity_backbone_kick_structure_validated|_gravity_backbone_kick_structure_validated]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:160-160`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:330-330`

**Downstream**

- `callees` → [[dynamics.perturbations_gravity_backbone_kick_acceleration_ii|gravity_backbone_kick_acceleration_ii]] · `callers` · call · `src/core/types/effector_sampling.jl:221-221`
<!-- vulcan:connections:end -->

## Limitations
Kick-based splitting is first-order accurate in the perturbation strength, and the hook offers no way to declare a preferred substep or ordering relative to other kicks. Torques are not supported by this path, so an effector with rotational effects must still use `wrench`. As with the other declaration hooks, the returned symbol is unchecked and a misspelled value falls through as unsupported without warning.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 218.
