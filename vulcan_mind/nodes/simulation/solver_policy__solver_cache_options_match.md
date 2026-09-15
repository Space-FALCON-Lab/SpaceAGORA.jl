---
id: simulation.solver_policy__solver_cache_options_match
label: _solver_cache_options_match
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solver_cache_options_match
  lines:
  - 271
  - 271
inputs:
- id: solver_cache
  type: SolverIntegratorCache
  units: n/a
  required: true
  description: Positional argument `solver_cache`.
- id: save_everystep
  type: Bool
  units: n/a
  required: true
  description: Positional argument `save_everystep`.
- id: save_on
  type: Bool
  units: n/a
  required: true
  description: Positional argument `save_on`.
- id: save_start
  type: Bool
  units: n/a
  required: true
  description: Positional argument `save_start`.
- id: save_end
  type: Bool
  units: n/a
  required: true
  description: Positional argument `save_end`.
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
  description: Return value of `_solver_cache_options_match`.
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

# _solver_cache_options_match

## Purpose
Checks that a cached integrator was initialised with exactly the save options the current solve resolved, gating safe reuse.

## Design & Implementation
Returns the conjunction of four `==` comparisons between `solver_cache.save_everystep/save_on/save_start/save_end` and the corresponding arguments. Called in `_solve_with_explicit_solver` before `reinit!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `solver_cache` | SolverIntegratorCache | n/a | yes | Positional argument `solver_cache`. |
| in | `save_everystep` | Bool | n/a | yes | Positional argument `save_everystep`. |
| in | `save_on` | Bool | n/a | yes | Positional argument `save_on`. |
| in | `save_start` | Bool | n/a | yes | Positional argument `save_start`. |
| in | `save_end` | Bool | n/a | yes | Positional argument `save_end`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_solver_cache_options_match`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:337-337`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only save flags are compared; algorithm, tolerances, `dtmax`, `maxiters`, and problem structure are assumed identical. A `nothing` integrator with matching flags still passes this check, so the caller must test `solver_cache.integrator !== nothing` separately, which it does.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 271.
