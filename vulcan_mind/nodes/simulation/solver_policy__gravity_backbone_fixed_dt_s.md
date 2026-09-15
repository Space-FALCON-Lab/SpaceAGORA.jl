---
id: simulation.solver_policy__gravity_backbone_fixed_dt_s
label: _gravity_backbone_fixed_dt_s
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _gravity_backbone_fixed_dt_s
  lines:
  - 93
  - 93
inputs:
- id: cfg
  type: SolverConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: Float64
  units: n/a
  description: Return value of `_gravity_backbone_fixed_dt_s`.
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

# _gravity_backbone_fixed_dt_s

## Purpose
Returns the macro step (seconds) for the gravity-backbone split integrator, defaulting to `dt_max_orbit` when `SolverConfig.gravity_backbone_dt_s` is unset.

## Design & Implementation
Mirrors `_symplectic_fixed_dt_s`: selects `cfg.gravity_backbone_dt_s` or `args.integration_tolerances.dt_max_orbit`, throws `ArgumentError` unless `dt > 0.0`, and offers a one-argument overload using `_active_solver_config()`. Consumed once per run by `_solve_with_gravity_backbone_solver`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_gravity_backbone_fixed_dt_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:565-565`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the same step is used for both the KahanLi8 core and the half-kick perturbation splitting, the perturbation accuracy is bounded by this single scalar with no adaptive refinement. `Inf` is accepted. The `dt_max_orbit` default is shared with unrelated adaptive modes.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 93.
