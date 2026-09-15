---
id: simulation.setup__ephemeris_explicit_cache_store_bang
label: _ephemeris_explicit_cache_store!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ephemeris_explicit_cache_store!
  lines:
  - 324
  - 324
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
- id: replace
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `replace` (default `false`).
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
  description: 'Return value of `_ephemeris_explicit_cache_store!`; mutates `cache`
    in place. Type parameters: `{K, T}`.'
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

# _ephemeris_explicit_cache_store!

## Purpose
Stores a cache in a reuse dictionary without an entry cap, used for explicitly registered (prewarmed) N-body caches where the operator has requested retention and optionally replacement.

## Design & Implementation
`_ephemeris_explicit_cache_store!(cache::AbstractDict{K,T}, key::K, value::T; replace::Bool=false)::T`. Inside `_EPHEMERIS_REUSE_LOCK`, when `replace` is false it returns any existing entry unchanged; when `replace` is true or no entry exists it assigns `cache[key] = value` and returns `value`. Called by `_register_prewarmed_nbody_ephemeris_cache!` against `_NBODY_EPHEMERIS_PREWARMED_CACHE`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | AbstractDict{K, T} | n/a | yes | Positional argument `cache`. |
| in | `key` | K | n/a | yes | Positional argument `key`. |
| in | `value` | T | n/a | yes | Positional argument `value`. |
| in | `replace` | Bool | n/a | no | Keyword argument `replace` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | T | n/a | — | Return value of `_ephemeris_explicit_cache_store!`; mutates `cache` in place. Type parameters: `{K, T}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__register_prewarmed_nbody_ephemeris_cache_bang|_register_prewarmed_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1590-1590`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No size bound at all, so repeated registrations with distinct keys grow the dictionary until `_clear_ephemeris_reuse_cache!` is called. Replacing an entry while another run holds a reference to the old cache is safe only because caches are immutable after build; nothing enforces that.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 324.
