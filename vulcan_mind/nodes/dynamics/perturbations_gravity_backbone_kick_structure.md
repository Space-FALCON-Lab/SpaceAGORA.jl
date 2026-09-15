---
id: dynamics.perturbations_gravity_backbone_kick_structure
label: gravity_backbone_kick_structure
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: gravity_backbone_kick_structure
  lines:
  - 974
  - 974
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
  description: Return value of `gravity_backbone_kick_structure`. Returns `:velocity_kick_explicit`.
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

# gravity_backbone_kick_structure

## Purpose
Tells the split solver that the N-body effector's contribution is applied as an explicit velocity kick rather than folded into the implicit backbone.

## Design & Implementation
Returns the symbol `:velocity_kick_explicit`. Declared `@inline`. The split partition queries this once at setup to decide which effectors go into the implicit static-gravity step and which are applied explicitly between steps.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `gravity_backbone_kick_structure`. Returns `:velocity_kick_explicit`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:209-209`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.solver_policy__gravity_backbone_kick_structure_validated|_gravity_backbone_kick_structure_validated]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:160-160`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:330-330`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Third-body forces vary slowly, so an explicit kick is accurate, but the classification is fixed and cannot be overridden per run.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 974.
