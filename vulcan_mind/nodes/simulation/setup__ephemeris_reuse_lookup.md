---
id: simulation.setup__ephemeris_reuse_lookup
label: _ephemeris_reuse_lookup
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ephemeris_reuse_lookup
  lines:
  - 300
  - 300
inputs:
- id: cache
  type: AbstractDict{K, T}
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: key
  type: K
  units: n/a
  required: true
  description: Positional argument `key`.
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
  description: 'Return value of `_ephemeris_reuse_lookup`. Returns `lock(_EPHEMERIS_REUSE_LOCK)
    do`. Type parameters: `{K, T}`.'
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

# _ephemeris_reuse_lookup

## Purpose
Thread-safe read from one of the process-global ephemeris reuse dictionaries, returning the cached object or `nothing`.

## Design & Implementation
Generic over `cache::AbstractDict{K, T}` and `key::K`. Acquires `_EPHEMERIS_REUSE_LOCK` (a `ReentrantLock`) with the `lock(f, l)` do-block form and evaluates `get(cache, key, nothing)` inside the critical section. Returns `Union{T, Nothing}`. Because the lock is reentrant, a caller already holding it (for example `_ephemeris_reuse_store!` calling back) does not deadlock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | AbstractDict{K, T} | n/a | yes | Positional argument `cache`. |
| in | `key` | K | n/a | yes | Positional argument `key`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_ephemeris_reuse_lookup`. Returns `lock(_EPHEMERIS_REUSE_LOCK) do`. Type parameters: `{K, T}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1827-1827`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1873-1873`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1769-1769`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1717-1717`
- [[simulation.setup__prewarmed_nbody_ephemeris_lookup|_prewarmed_nbody_ephemeris_lookup]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1577-1577`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Returning the cached object itself, not a copy, means callers share mutable cache state across runs; the caches must therefore be treated as read-only after construction. The lock is global across all three cache kinds, so concurrent setups of unrelated caches serialise on it.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 300.
