---
id: simulation.setup__ephemeris_reuse_enabled
label: _ephemeris_reuse_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _ephemeris_reuse_enabled
  lines:
  - 232
  - 232
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
  type: Bool
  units: n/a
  description: Return value of `_ephemeris_reuse_enabled`.
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

# _ephemeris_reuse_enabled

## Purpose
Controls whether ephemeris caches built for one run are retained in process-global dictionaries and handed to later runs with an identical key, avoiding repeated SPICE sampling across a sweep.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EPHEMERIS_CACHE_REUSE", true)`. When true, the `_initialize_*_ephemeris_cache!` functions call `_ephemeris_reuse_lookup` before building and `_ephemeris_reuse_store!` after, using `_SRP_EPHEMERIS_REUSE_CACHE`, `_NBODY_EPHEMERIS_REUSE_CACHE`, and `_PLANET_FRAME_EPHEMERIS_REUSE_CACHE`. Throws `ArgumentError` on malformed values.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_ephemeris_reuse_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1825-1825`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1871-1871`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1767-1767`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:233-233`
<!-- vulcan:connections:end -->

## Limitations
Reuse keys quantise times to microseconds and planet constants to fixed scales, so two runs differing below those resolutions share a cache. The global dictionaries are never cleared automatically except via `_clear_ephemeris_reuse_cache!`, so a long-lived process accumulates up to `_ephemeris_reuse_max_entries` caches per kind.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 232.
