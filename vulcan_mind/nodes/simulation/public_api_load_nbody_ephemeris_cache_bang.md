---
id: simulation.public_api_load_nbody_ephemeris_cache_bang
label: load_nbody_ephemeris_cache!
kind: function
source:
  file: src/simulation/engine/public_api.jl
  symbol: load_nbody_ephemeris_cache!
  lines:
  - 53
  - 53
inputs:
- id: path
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `path`.
- id: replace
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `replace` (default `true`).
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
  description: Return value of `load_nbody_ephemeris_cache!`; mutates `path` in place.
    Returns `_load_nbody_ephemeris_cache!(String(path); replace=replace)`.
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

# load_nbody_ephemeris_cache!

## Purpose
Loads a previously saved n-body ephemeris cache from disk into the process-global cache, letting repeated simulations skip the prewarm step entirely.

## Design & Implementation
The signature is `load_nbody_ephemeris_cache!(path::AbstractString; replace::Bool=true)`. It normalises the argument with `String(path)` so that substrings and other `AbstractString` views do not leak into the cache implementation, then delegates to `_load_nbody_ephemeris_cache!`. With `replace=true`, the default, the existing cache contents are discarded in favour of the file; with `replace=false` the loaded entries are merged into whatever is already resident.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | AbstractString | n/a | yes | Positional argument `path`. |
| in | `replace` | Bool | n/a | no | Keyword argument `replace` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `load_nbody_ephemeris_cache!`; mutates `path` in place. Returns `_load_nbody_ephemeris_cache!(String(path); replace=replace)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/public_api.jl`
- [[spaceagora.spaceagora_prewarm_nbody_ephemeris_cache|prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:536-536`

**Downstream**

- `callees` → [[simulation.setup__load_nbody_ephemeris_cache_bang|_load_nbody_ephemeris_cache!]] · `callers` · call · `src/simulation/engine/public_api.jl:54-54`
<!-- vulcan:connections:end -->

## Limitations
The trailing bang refers to mutation of shared, process-wide state, so concurrent calls from several tasks race against each other and against any in-flight integration that is reading the cache. Nothing at this layer checks that the cached epochs, step size, or body set match the simulation about to run, so a stale or mismatched file can silently supply the wrong third-body positions. Missing or malformed files surface as errors from the inner loader.

## Provenance
Mapped from `src/simulation/engine/public_api.jl` line 53.
