---
id: simulation.setup__ephemeris_reuse_store_bang
label: _ephemeris_reuse_store!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ephemeris_reuse_store!
  lines:
  - 306
  - 306
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
- id: value
  type: T
  units: n/a
  required: true
  description: Positional argument `value`.
- id: max_entries
  type: Int
  units: n/a
  required: true
  description: Positional argument `max_entries`.
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
  type: T
  units: n/a
  description: 'Return value of `_ephemeris_reuse_store!`; mutates `cache` in place.
    Type parameters: `{K, T}`.'
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

# _ephemeris_reuse_store!

## Purpose
Inserts a freshly built ephemeris cache into a reuse dictionary unless an equal key is already present, enforcing the entry cap by evicting one existing entry.

## Design & Implementation
`_ephemeris_reuse_store!(cache::AbstractDict{K,T}, key::K, value::T, max_entries::Int)::T`. Under `_EPHEMERIS_REUSE_LOCK` it first returns any `existing` entry for `key` (first-writer wins, so concurrent builders converge on one object). If `max_entries <= 0` the value is returned without insertion. If `length(cache) >= max_entries` it deletes `first(keys(cache))` before assigning `cache[key] = value`. Always returns the object that is now canonical for the key. Mutates `cache`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | AbstractDict{K, T} | n/a | yes | Positional argument `cache`. |
| in | `key` | K | n/a | yes | Positional argument `key`. |
| in | `value` | T | n/a | yes | Positional argument `value`. |
| in | `max_entries` | Int | n/a | yes | Positional argument `max_entries`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | T | n/a | — | Return value of `_ephemeris_reuse_store!`; mutates `cache` in place. Type parameters: `{K, T}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1844-1844`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1895-1895`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1790-1790`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`first(keys(cache))` on a `Dict` is hash-order, not insertion order, so eviction is arbitrary rather than LRU or FIFO. Only one entry is evicted per call, which is sufficient only because insertions are one at a time. When an existing entry is returned the caller's freshly built `value` is discarded, wasting the sampling work that produced it.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 306.
