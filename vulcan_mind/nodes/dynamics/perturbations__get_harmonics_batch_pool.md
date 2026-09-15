---
id: dynamics.perturbations__get_harmonics_batch_pool
label: _get_harmonics_batch_pool
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _get_harmonics_batch_pool
  lines:
  - 335
  - 335
inputs:
- id: model
  type: GravitationalHarmonicsModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: n_workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_workers`.
- id: batch_size
  type: Int
  units: n/a
  required: true
  description: Positional argument `batch_size`.
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
  type: Vector{HarmonicsBatchWorkspace}
  units: n/a
  description: Return value of `_get_harmonics_batch_pool`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _get_harmonics_batch_pool

## Purpose
Returns a per-worker pool of harmonics batch workspaces for a model, creating or enlarging it under a lock when the current pool is too small.

## Design & Implementation
Keys `_HARMONICS_BATCH_POOL` by `objectid(model)`. A fast unlocked read returns the pool if it has at least `n_workers` entries sized to at least `batch_size`; otherwise the lock is taken, the check repeated, and a fresh vector of workspaces built. Returns the pool.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `n_workers` | Int | n/a | yes | Positional argument `n_workers`. |
| in | `batch_size` | Int | n/a | yes | Positional argument `batch_size`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{HarmonicsBatchWorkspace} | n/a | — | Return value of `_get_harmonics_batch_pool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__get_harmonics_batch_pool_cached_bang|_get_harmonics_batch_pool_cached!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:368-368`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__make_harmonics_batch_workspace|_make_harmonics_batch_workspace]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:348-348`
<!-- vulcan:connections:end -->

## Limitations
Keyed by object identity, so a model garbage-collected and reallocated at the same address would alias a stale pool; and the process-global dictionary never evicts.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 335.
