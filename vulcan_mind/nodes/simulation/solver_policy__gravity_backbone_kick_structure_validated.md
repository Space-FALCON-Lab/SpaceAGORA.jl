---
id: simulation.solver_policy__gravity_backbone_kick_structure_validated
label: _gravity_backbone_kick_structure_validated
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _gravity_backbone_kick_structure_validated
  lines:
  - 159
  - 159
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
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
  description: Return value of `_gravity_backbone_kick_structure_validated`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _gravity_backbone_kick_structure_validated

## Purpose
Queries an effector's `gravity_backbone_kick_structure` and enforces that it is `:unsupported` or `:velocity_kick_explicit`.

## Design & Implementation
Calls `SimulationModel.gravity_backbone_kick_structure(effector)`; returns the symbol when it matches one of the two allowed values, otherwise throws `ArgumentError` with the effector type name and `repr(structure)`. Used by `_gravity_backbone_has_kicks` and `_gravity_backbone_reject_reason`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_gravity_backbone_kick_structure_validated`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.dynamics_rhs__gravity_backbone_kick_acceleration|_gravity_backbone_kick_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1516-1516`
- [[simulation.solver_policy__gravity_backbone_has_kicks|_gravity_backbone_has_kicks]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:171-171`
- [[simulation.solver_policy__gravity_backbone_reject_reason|_gravity_backbone_reject_reason]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:193-193`

**Downstream**

- `callees` → [[core.effector_sampling_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `callers` · call · `src/simulation/engine/solver_policy.jl:160-160`
- `callees` → [[dynamics.perturbations_gravity_backbone_kick_structure|gravity_backbone_kick_structure]] · `callers` · call · `src/simulation/engine/solver_policy.jl:160-160`
<!-- vulcan:connections:end -->

## Limitations
Only explicit velocity kicks are supported; position-dependent or implicit perturbations cannot be expressed. Errors are raised eagerly during eligibility scanning rather than deferred to solve time.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 159.
