---
id: simulation.vacuum_predicted_gram__vacuum_gram_cache_for_sat_bang
label: _vacuum_gram_cache_for_sat!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _vacuum_gram_cache_for_sat!
  lines:
  - 46
  - 46
inputs:
- id: caches
  type: Vector{Union{Nothing, VacuumPredictedGRAMCache}}
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
  type: VacuumPredictedGRAMCache
  units: n/a
  description: Return value of `_vacuum_gram_cache_for_sat!`; mutates `caches` in
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

# _vacuum_gram_cache_for_sat!

## Purpose
Fetches, and lazily creates, the `VacuumPredictedGRAMCache` belonging to satellite `sat_idx` from the per-run vector of caches, so that each spacecraft in a constellation keeps an independent vacuum prediction and spline set.

## Design & Implementation
Signature `_vacuum_gram_cache_for_sat!(caches::Vector{Union{Nothing, VacuumPredictedGRAMCache}}, sat_idx::Int)::VacuumPredictedGRAMCache`, marked `@inline`. If `sat_idx <= length(caches)` it reads the slot with `@inbounds`; when the slot holds `nothing` it constructs a fresh `VacuumPredictedGRAMCache()` and stores it back into `caches[sat_idx]` (the mutation that earns the `!`). The cache object is returned. If `sat_idx` exceeds the vector length, a new temporary cache is returned without being stored, so the caller still gets a usable object.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `caches` | Vector{Union{Nothing, VacuumPredictedGRAMCache}} | n/a | yes | Positional argument `caches`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VacuumPredictedGRAMCache | n/a | — | Return value of `_vacuum_gram_cache_for_sat!`; mutates `caches` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:102-102`

**Downstream**

- `callees` → [[core.runtime_types_vacuumpredictedgramcache|VacuumPredictedGRAMCache]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:53-53`
- `callees` → [[simulation.vacuum_predicted_gram_vacuumpredictedgramcache|VacuumPredictedGRAMCache]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:53-53`
<!-- vulcan:connections:end -->

## Limitations
The out-of-range branch allocates a throwaway cache on every call, which defeats caching entirely and silently hides a sizing bug in the caller that allocated `caches` too short. A `sat_idx` of zero or negative passes the length test and then triggers an out-of-bounds read that `@inbounds` turns into undefined behaviour instead of a `BoundsError`. Two threads creating the cache for the same satellite simultaneously can race and one cache is lost.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl` line 46.
