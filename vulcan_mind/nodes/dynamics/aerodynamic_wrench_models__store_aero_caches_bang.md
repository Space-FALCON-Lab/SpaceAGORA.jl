---
id: dynamics.aerodynamic_wrench_models__store_aero_caches_bang
label: _store_aero_caches!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _store_aero_caches!
  lines:
  - 198
  - 198
inputs:
- id: param
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `param`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: drag_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `drag_ii`.
- id: lift_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `lift_ii`.
- id: cross_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `cross_ii`.
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
  type: Nothing
  units: n/a
  description: Return value of `_store_aero_caches!`; mutates `param` in place.
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

# _store_aero_caches!

## Purpose
Records per-satellite inertial drag, lift, and cross-force vectors into `param.save_cache` so output callbacks can log them without recomputing aerodynamics.

## Design & Implementation
Three calls to `_store_vector_cache!` targeting `param.save_cache.drag_cache`, `lift_cache`, and `cross_cache` with the respective `SVector{3,Float64}` inputs at index `sat_idx`. Returns `nothing`. Invoked by every `wrench_caching!` method and by the fM `calcForceTorque`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `drag_ii` | SVector{3, Float64} | n/a | yes | Positional argument `drag_ii`. |
| in | `lift_ii` | SVector{3, Float64} | n/a | yes | Positional argument `lift_ii`. |
| in | `cross_ii` | SVector{3, Float64} | n/a | yes | Positional argument `cross_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_store_aero_caches!`; mutates `param` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:721-721`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:890-890`
- [[dynamics.aerodynamic_wrench_models_wrench_caching_bang|wrench_caching!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:553-553`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models__store_vector_cache_bang|_store_vector_cache!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:205-205`
<!-- vulcan:connections:end -->

## Limitations
Caches store only the most recent evaluation, so with adaptive solvers that reject steps the cached values may correspond to a rejected trial state, not the accepted step. Not thread-safe across satellites if caches are grown concurrently.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 198.
