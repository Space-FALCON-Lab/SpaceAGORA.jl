---
id: simulation.model_selection__gram_density_cache_for_sat_bang
label: _gram_density_cache_for_sat!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _gram_density_cache_for_sat!
  lines:
  - 115
  - 115
inputs:
- id: caches
  type: Vector{Union{Nothing, GramTrackCache}}
  units: n/a
  required: true
  description: Positional argument `caches`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: GramTrackCache
  units: n/a
  description: Return value of `_gram_density_cache_for_sat!`; mutates `caches` in
    place.
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

# _gram_density_cache_for_sat!

## Purpose
Returns the `GramTrackCache` associated with satellite `sat_idx`, lazily allocating and storing one in the shared `caches` vector on first use so that GRAM track lookups for each satellite persist across integrator steps.

## Design & Implementation
Takes `caches::Vector{Union{Nothing, GramTrackCache}}` and `sat_idx::Int`. When `sat_idx <= length(caches)`, it reads `caches[sat_idx]` with `@inbounds`; if the slot is `nothing` a fresh `GramTrackCache()` is constructed and written back into the slot (this is the mutation implied by the `!` suffix). When `sat_idx` exceeds the vector length a throwaway `GramTrackCache()` is returned without storing it. The return annotation `::GramTrackCache` guarantees a non-`nothing` result to the caller.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `caches` | Vector{Union{Nothing, GramTrackCache}} | n/a | yes | Positional argument `caches`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GramTrackCache | n/a | — | Return value of `_gram_density_cache_for_sat!`; mutates `caches` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:125-125`

**Downstream**

- `callees` → [[core.runtime_types_gramtrackcache|GramTrackCache]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:122-122`
- `callees` → [[simulation.config_gramtrackcache|GramTrackCache]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:122-122`
<!-- vulcan:connections:end -->

## Limitations
Out-of-range satellites get a new cache every call, silently defeating caching and allocating on each density evaluation. Negative or zero `sat_idx` passes the upper-bound check and triggers an `@inbounds` read at an invalid index, which is undefined behaviour rather than a clean error. Lazy initialisation is not thread-safe: two workers evaluating the same `sat_idx` concurrently could both allocate and one write would be lost.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl` line 115.
