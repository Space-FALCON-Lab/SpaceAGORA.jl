---
id: simulation.solver_policy_solverintegratorcache
label: SolverIntegratorCache
kind: struct
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: SolverIntegratorCache
  lines:
  - 257
  - 257
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Field `integrator`.
- id: save_everystep
  type: Bool
  units: n/a
  required: true
  description: Field `save_everystep`.
- id: save_on
  type: Bool
  units: n/a
  required: true
  description: Field `save_on`.
- id: save_start
  type: Bool
  units: n/a
  required: true
  description: Field `save_start`.
- id: save_end
  type: Bool
  units: n/a
  required: true
  description: Field `save_end`.
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
  type: SolverIntegratorCache
  units: n/a
  description: Constructed `SolverIntegratorCache`.
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

# SolverIntegratorCache

## Purpose
Mutable holder for a reusable DiffEq integrator plus the four save options it was initialised with, so repeated runs can `reinit!` instead of re-allocating.

## Design & Implementation
Fields: `integrator::Any` (starts `nothing`), and `save_everystep`, `save_on`, `save_start`, `save_end` (all `Bool`, initialised `true`). The default constructor `SolverIntegratorCache()` sets those values. `_solver_cache_options_match` compares a request against the stored flags and `_cache_integrator!` overwrites all five fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Field `integrator`. |
| in | `save_everystep` | Bool | n/a | yes | Field `save_everystep`. |
| in | `save_on` | Bool | n/a | yes | Field `save_on`. |
| in | `save_start` | Bool | n/a | yes | Field `save_start`. |
| in | `save_end` | Bool | n/a | yes | Field `save_end`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SolverIntegratorCache | n/a | — | Constructed `SolverIntegratorCache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`integrator::Any` makes every access type-unstable. The cache does not record the algorithm, tolerances, `dtmax`, or `maxiters`, so reusing it with different values silently keeps the original settings. Not thread-safe; concurrent runs must own separate caches.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 257.
