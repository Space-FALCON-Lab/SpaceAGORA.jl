---
id: dynamics.aerodynamic_wrench_models__store_vector_cache_bang
label: _store_vector_cache!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _store_vector_cache!
  lines:
  - 181
  - 181
inputs:
- id: cache
  type: Vector{SVector{3, Float64}}
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: value
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `value`.
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
  description: Return value of `_store_vector_cache!`; mutates `cache` in place.
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

# _store_vector_cache!

## Purpose
Writes a 3-vector into slot `sat_idx` of a save cache, growing the vector and zero-initialising any unassigned entries first.

## Design & Implementation
If `length(cache) < sat_idx` it `resize!`s to `sat_idx` and loops all indices, assigning zero `SVector` where `!isassigned`. Then `cache[sat_idx] = value` under `@inbounds`. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | Vector{SVector{3, Float64}} | n/a | yes | Positional argument `cache`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `value` | SVector{3, Float64} | n/a | yes | Positional argument `value`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_store_vector_cache!`; mutates `cache` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__store_aero_caches_bang|_store_aero_caches!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:205-205`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `SVector{3,Float64}` is `isbits`, `resize!` leaves new slots holding garbage bits that `isassigned` reports as assigned, so the zero-fill loop is ineffective for those entries; unwritten satellites may expose uninitialised values until first written. The resize is a growth-only path and never shrinks.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 181.
