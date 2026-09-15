---
id: simulation.setup__ephemeris_reuse_max_entries
label: _ephemeris_reuse_max_entries
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ephemeris_reuse_max_entries
  lines:
  - 236
  - 236
inputs:
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
  type: Int
  units: n/a
  description: Return value of `_ephemeris_reuse_max_entries`.
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

# _ephemeris_reuse_max_entries

## Purpose
Bounds how many distinct ephemeris caches of each kind the process-global reuse dictionaries retain, preventing unbounded memory growth in long parameter sweeps.

## Design & Implementation
Returns `_parse_nonnegative_int_env("SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES", 32)`. Zero is legal and makes `_ephemeris_reuse_store!` skip insertion entirely (the built cache is still returned to the caller). Negative or non-integer text throws `ArgumentError`. The same limit applies independently to the SRP, N-body and planet-frame dictionaries.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_ephemeris_reuse_max_entries`. |
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

- `callees` → [[simulation.setup__parse_nonnegative_int_env|_parse_nonnegative_int_env]] · `callers` · call · `src/simulation/engine/setup.jl:237-237`
<!-- vulcan:connections:end -->

## Limitations
Eviction in `_ephemeris_reuse_store!` removes `first(keys(cache))`, which for a `Dict` is not insertion order, so the policy is effectively arbitrary rather than LRU. The limit is a count, not a byte budget, so 32 large N-body caches can still consume gigabytes.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 236.
