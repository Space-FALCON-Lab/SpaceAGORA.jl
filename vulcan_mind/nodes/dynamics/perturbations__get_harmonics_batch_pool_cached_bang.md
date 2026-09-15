---
id: dynamics.perturbations__get_harmonics_batch_pool_cached_bang
label: _get_harmonics_batch_pool_cached!
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _get_harmonics_batch_pool_cached!
  lines:
  - 356
  - 356
inputs:
- id: pool_ref
  type: Base.RefValue{Any}
  units: n/a
  required: true
  description: Positional argument `pool_ref`.
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
  description: Return value of `_get_harmonics_batch_pool_cached!`; mutates `pool_ref`
    in place.
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

# _get_harmonics_batch_pool_cached!

## Purpose
The per-run fast path for the batch pool: reads a `Ref` on shared buffers before falling back to the global dictionary.

## Design & Implementation
If `pool_ref[]` is a workspace vector large enough for `n_workers` and `batch_size`, returns it; otherwise fetches from `_get_harmonics_batch_pool`, stores it in the ref, and returns it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pool_ref` | Base.RefValue{Any} | n/a | yes | Positional argument `pool_ref`. |
| in | `model` | GravitationalHarmonicsModel | n/a | yes | Positional argument `model`. |
| in | `n_workers` | Int | n/a | yes | Positional argument `n_workers`. |
| in | `batch_size` | Int | n/a | yes | Positional argument `batch_size`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{HarmonicsBatchWorkspace} | n/a | — | Return value of `_get_harmonics_batch_pool_cached!`; mutates `pool_ref` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`
- [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:944-944`

**Downstream**

- `callees` → [[dynamics.perturbations__get_harmonics_batch_pool|_get_harmonics_batch_pool]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:368-368`
<!-- vulcan:connections:end -->

## Limitations
The ref is `Ref{Any}`, so the `isa` check is a dynamic type test on every RHS call.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 356.
