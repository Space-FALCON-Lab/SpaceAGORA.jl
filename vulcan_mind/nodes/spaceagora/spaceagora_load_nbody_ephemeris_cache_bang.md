---
id: spaceagora.spaceagora_load_nbody_ephemeris_cache_bang
label: load_nbody_ephemeris_cache!
kind: function
source:
  file: src/SpaceAGORA.jl
  symbol: load_nbody_ephemeris_cache!
  lines:
  - 544
  - 544
inputs:
- id: args
  type: Vararg{Any}
  units: n/a
  required: false
  description: Positional argument `args` (variadic).
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
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
  description: Return value of `load_nbody_ephemeris_cache!`; mutates `args` in place.
    Returns `SimulationEngine.load_nbody_ephemeris_cache!(args...; kwargs...)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- spaceagora
charts:
- spaceagora
origin: agent
---

# load_nbody_ephemeris_cache!

## Purpose
`load_nbody_ephemeris_cache!` deserializes an N-body ephemeris cache previously written by `prewarm_nbody_ephemeris_cache(...; save_path=...)` and registers it in the current Julia process, so worker processes in a multi-process Monte Carlo campaign can share one precomputed SPICE sample set. The package-level definition forwards to `SimulationEngine.load_nbody_ephemeris_cache!`.

## Design & Implementation
The concrete method is `load_nbody_ephemeris_cache!(path::AbstractString; replace::Bool=true)`, which converts `path` to `String` and delegates to `_load_nbody_ephemeris_cache!(path; replace)`. The `!` marks that it mutates process-global state: the loaded cache is inserted into the engine's cache registry, overwriting any entry with the same deterministic key when `replace` is `true`. The cache is returned so callers can inspect its key or sample count. It is intended to be called once per worker before the first `run_simulation`, typically from a `@everywhere` block.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `load_nbody_ephemeris_cache!`; mutates `args` in place. Returns `SimulationEngine.load_nbody_ephemeris_cache!(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.spaceagora|SpaceAGORA]] · `api` → `module_api` · call · `src/SpaceAGORA.jl`
- [[spaceagora.spaceagora_prewarm_nbody_ephemeris_cache|prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:536-536`

**Downstream**

- `callees` → [[misc.assets_check_assets|check_assets]] · `callers` · call · `src/SpaceAGORA.jl:547-547`
- `callees` → [[spaceagora.spaceagora_check_assets|check_assets]] · `callers` · call · `src/SpaceAGORA.jl:547-547`
<!-- vulcan:connections:end -->

## Limitations
The file format is whatever the engine's serializer produced; loading a cache written by a different package version or Julia version can fail or produce a subtly incompatible object with no version check exposed here. A missing or unreadable `path` throws from the underlying I/O. With `replace=false` an existing entry is kept silently, so a stale cache can persist without warning. Registration is not thread-synchronised at this layer, so concurrent loads on the same process are a hazard.

## Provenance
Mapped from `src/SpaceAGORA.jl` line 544.
