---
id: simulation.solver_policy__cache_integrator_bang
label: _cache_integrator!
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _cache_integrator!
  lines:
  - 284
  - 284
inputs:
- id: solver_cache
  type: SolverIntegratorCache
  units: n/a
  required: true
  description: Positional argument `solver_cache`.
- id: integ
  type: Any
  units: n/a
  required: true
  description: Positional argument `integ`.
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
  type: Any
  units: n/a
  description: Return value of `_cache_integrator!`; mutates `solver_cache` in place.
    Returns `integ`.
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

# _cache_integrator!

## Purpose
Stores a freshly initialised integrator and its save options into a `SolverIntegratorCache`, returning the integrator for immediate use.

## Design & Implementation
Assigns `solver_cache.integrator = integ` and the four `Bool` save fields, then returns `integ`. Called from both the `maxiters === nothing` and `maxiters` branches of `_solve_with_explicit_solver`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `solver_cache` | SolverIntegratorCache | n/a | yes | Positional argument `solver_cache`. |
| in | `integ` | Any | n/a | yes | Positional argument `integ`. |
| in | `save_everystep` | Bool | n/a | yes | Positional argument `save_everystep`. |
| in | `save_on` | Bool | n/a | yes | Positional argument `save_on`. |
| in | `save_start` | Bool | n/a | yes | Positional argument `save_start`. |
| in | `save_end` | Bool | n/a | yes | Positional argument `save_end`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_cache_integrator!`; mutates `solver_cache` in place. Returns `integ`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:351-351`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Overwrites any previously cached integrator without releasing it explicitly, so two large integrators may coexist until garbage collection. Does not record the algorithm or tolerances, leaving mismatched reuse undetectable.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 284.
