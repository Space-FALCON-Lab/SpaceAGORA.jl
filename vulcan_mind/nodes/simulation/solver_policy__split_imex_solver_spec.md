---
id: simulation.solver_policy__split_imex_solver_spec
label: _split_imex_solver_spec
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _split_imex_solver_spec
  lines:
  - 215
  - 215
inputs:
- id: cfg
  type: SolverConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: Any
  units: n/a
  description: Return value of `_split_imex_solver_spec`.
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

# _split_imex_solver_spec

## Purpose
Maps `SolverConfig.split_imex_solver` to a concrete IMEX Runge-Kutta algorithm instance and label for `:split_imex` mode.

## Design & Implementation
Returns `(alg=KenCarp4(autodiff=AutoFiniteDiff()), label="KenCarp4")` for `:kencarp4`, and analogous tuples for `:kencarp47` and `:kencarp58`. Any other symbol throws `ArgumentError` listing the three supported values. A zero-argument overload reads `_active_solver_config()`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_split_imex_solver_spec`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__rodas5p_alg|_rodas5p_alg]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:710-710`
- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:710-710`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
All variants use finite-difference Jacobians, so the implicit stage cost scales with state dimension and no sparse `jac_prototype` is passed, unlike `_rodas5p_alg`. A new algorithm instance is constructed on every call.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 215.
