---
id: simulation.solver_policy__rodas5p_alg
label: _rodas5p_alg
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _rodas5p_alg
  lines:
  - 669
  - 669
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
  type: Any
  units: n/a
  description: Return value of `_rodas5p_alg`. Returns `prob.f.jac_prototype !== nothing
    ?`.
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

# _rodas5p_alg

## Purpose
Local closure inside `_solve_with_solver_policy` that constructs the `Rodas5P` implicit algorithm, choosing a sparse KLU linear solver only when the problem supplies a `jac_prototype`.

## Design & Implementation
Defined as `_rodas5p_alg() = prob.f.jac_prototype !== nothing ? Rodas5P(autodiff=AutoFiniteDiff(), linsolve=KLUFactorization()) : Rodas5P(autodiff=AutoFiniteDiff())`. It closes over `prob` and is used for the `:rodas5p` mode directly and as the stiff member of `AutoTsit5` in `:auto_stiff` mode. The comment explains that KLU on a dense W matrix gives wrong Newton corrections, hence the guard.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rodas5p_alg`. Returns `prob.f.jac_prototype !== nothing ?`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:669-669`

**Downstream**

- `callees` → [[simulation.solver_policy__auto_stiff_switched|_auto_stiff_switched]] · `callers` · call · `src/simulation/engine/solver_policy.jl:700-700`
- `callees` → [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:674-674`
- `callees` → [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:721-721`
- `callees` → [[simulation.solver_policy__split_imex_solver_spec|_split_imex_solver_spec]] · `callers` · call · `src/simulation/engine/solver_policy.jl:710-710`
- `callees` → [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callers` · feedback · `src/simulation/engine/solver_policy.jl:750-750`
<!-- vulcan:connections:end -->

## Limitations
Finite-difference Jacobians (`AutoFiniteDiff`) cost one RHS evaluation per state column per Jacobian, which dominates implicit-step cost for large multi-satellite states even with sparsity. The closure is redefined on every call to `_solve_with_solver_policy`. Presence of `jac_prototype` is treated as proof of sparsity without checking its density.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 669.
