---
id: simulation.solver_policy__gravity_backbone_has_kicks
label: _gravity_backbone_has_kicks
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _gravity_backbone_has_kicks
  lines:
  - 169
  - 169
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Bool
  units: n/a
  description: Return value of `_gravity_backbone_has_kicks`.
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

# _gravity_backbone_has_kicks

## Purpose
Reports whether any dynamic effector contributes an explicit velocity kick, which changes the gravity-backbone solver label and enables the half-kick stages.

## Design & Implementation
Iterates `args.dynamics_model.dynamic_effectors` with `@inbounds` and returns `true` at the first effector whose `_gravity_backbone_kick_structure_validated` result is `:velocity_kick_explicit`; otherwise `false`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gravity_backbone_has_kicks`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:552-552`

**Downstream**

- `callees` → [[simulation.solver_policy__gravity_backbone_kick_structure_validated|_gravity_backbone_kick_structure_validated]] · `callers` · call · `src/simulation/engine/solver_policy.jl:171-171`
<!-- vulcan:connections:end -->

## Limitations
Validation throws for effectors with malformed kick structure even though this function only needs a boolean. The result is only used for the label string in `_solve_with_gravity_backbone_solver`; the half-kick calls run unconditionally whenever `half_dt > 0`, relying on `_gravity_backbone_half_kick!` to be a no-op without kicks.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 169.
