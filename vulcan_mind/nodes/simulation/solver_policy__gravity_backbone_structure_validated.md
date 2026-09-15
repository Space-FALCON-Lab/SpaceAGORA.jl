---
id: simulation.solver_policy__gravity_backbone_structure_validated
label: _gravity_backbone_structure_validated
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _gravity_backbone_structure_validated
  lines:
  - 149
  - 149
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
  description: Return value of `_gravity_backbone_structure_validated`.
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

# _gravity_backbone_structure_validated

## Purpose
Queries an effector's `gravity_backbone_structure` classification and enforces that it is one of the two recognised values.

## Design & Implementation
Calls `SimulationModel.gravity_backbone_structure(effector)` and returns the symbol if it is `:unsupported` or `:position_only_static_gravity`; any other value throws `ArgumentError` naming the effector type and the offending symbol via `repr`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_gravity_backbone_structure_validated`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.dynamics_rhs__gravity_backbone_core_acceleration|_gravity_backbone_core_acceleration]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1499-1499`
- [[simulation.solver_policy__gravity_backbone_reject_reason|_gravity_backbone_reject_reason]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:192-192`

**Downstream**

- `callees` → [[core.effector_sampling_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/simulation/engine/solver_policy.jl:150-150`
- `callees` → [[dynamics.perturbations_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/simulation/engine/solver_policy.jl:150-150`
- `callees` → [[environment.gravity_models_gravity_backbone_structure|gravity_backbone_structure]] · `callers` · call · `src/simulation/engine/solver_policy.jl:150-150`
<!-- vulcan:connections:end -->

## Limitations
A misbehaving effector method turns a configuration query into a hard error even when the backbone mode is not selected; `_gravity_backbone_reject_reason` calls this for every effector, so one bad implementation blocks the whole eligibility check. The accepted vocabulary is duplicated as literals here and in the kick validator.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 149.
