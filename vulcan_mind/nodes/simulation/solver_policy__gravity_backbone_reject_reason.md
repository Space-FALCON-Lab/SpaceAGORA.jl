---
id: simulation.solver_policy__gravity_backbone_reject_reason
label: _gravity_backbone_reject_reason
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _gravity_backbone_reject_reason
  lines:
  - 183
  - 183
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_gravity_backbone_reject_reason`.
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

# _gravity_backbone_reject_reason

## Purpose
Returns a descriptive reason the `gravity_backbone_split` solver mode cannot be used for this run, or `nothing` when every effector fits the core-plus-kick structure.

## Design & Implementation
Rejects on `orientation_sim`, any control/guidance/navigation effectors, or no dynamic effectors. For each effector it validates both `core_structure` and `kick_structure`; an effector must be either `:position_only_static_gravity` (sets `has_core`) or `:velocity_kick_explicit`, else it is rejected. Core effectors may not require atmosphere, solar, or third-body ephemeris (`req.third_body_names` non-empty); kick effectors may not require atmosphere or `planet_frame`. Finally requires `has_core`. All messages are prefixed with `SPACEAGORA_SOLVER_MODE=gravity_backbone_split`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_gravity_backbone_reject_reason`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:652-652`

**Downstream**

- `callees` → [[core.effector_sampling_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:199-199`
- `callees` → [[dynamics.aerodynamic_wrench_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:199-199`
- `callees` → [[dynamics.perturbations_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:199-199`
- `callees` → [[environment.gravity_models_environment_requirements|environment_requirements]] · `callers` · call · `src/simulation/engine/solver_policy.jl:199-199`
- `callees` → [[simulation.solver_policy__gravity_backbone_kick_structure_validated|_gravity_backbone_kick_structure_validated]] · `callers` · call · `src/simulation/engine/solver_policy.jl:193-193`
- `callees` → [[simulation.solver_policy__gravity_backbone_structure_validated|_gravity_backbone_structure_validated]] · `callers` · call · `src/simulation/engine/solver_policy.jl:192-192`
<!-- vulcan:connections:end -->

## Limitations
The messages reference the environment variable name even when the mode came from a `SolverConfig` struct. Structure validators throw rather than return a reason for malformed effectors. Kick effectors that need solar samples are allowed, so SRP kicks are permitted while a solar-dependent gravity core is not, which is asymmetric but intentional.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 183.
