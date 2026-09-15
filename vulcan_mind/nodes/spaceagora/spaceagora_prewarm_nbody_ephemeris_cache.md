---
id: spaceagora.spaceagora_prewarm_nbody_ephemeris_cache
label: prewarm_nbody_ephemeris_cache
kind: function
source:
  file: src/SpaceAGORA.jl
  symbol: prewarm_nbody_ephemeris_cache
  lines:
  - 533
  - 533
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
  type: SimulationEngine.prewarm_nbody_ephemeris_cache
  units: n/a
  description: Return value of `prewarm_nbody_ephemeris_cache`. Returns `SimulationEngine.prewarm_nbody_ephemeris_cache(args...;
    kwargs...)`.
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

# prewarm_nbody_ephemeris_cache

## Purpose
`prewarm_nbody_ephemeris_cache` precomputes and registers a process-local N-body SPICE ephemeris cache so that subsequent `run_simulation` calls in a Monte Carlo campaign can reuse the same sampled third-body positions instead of querying SPICE per run. The package-level symbol forwards to `SimulationEngine.prewarm_nbody_ephemeris_cache`.

## Design & Implementation
Two methods exist in `SimulationEngine`: `prewarm_nbody_ephemeris_cache(args; dt_s::Union{Nothing,Real}=nothing, mission_end_s::Union{Nothing,Real}=nothing, save_path::Union{Nothing,AbstractString}=nothing)` calls the private `_prewarm_nbody_ephemeris_cache`, while `prewarm_nbody_ephemeris_cache(config::SimulationEngineConfig, args; kwargs...)` wraps the same call in `_with_engine_env_overrides(config, ...)` so the config's environment settings apply during the build. `dt_s` is the cache sample spacing in seconds and `mission_end_s` the span to cover; when `nothing`, the runtime derives them from `args` using the same deterministic key the simulation setup uses. If `save_path` is given the cache is also serialized to disk for `load_nbody_ephemeris_cache!` on worker processes. The cache object is returned and registered in the current process.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationEngine.prewarm_nbody_ephemeris_cache | n/a | — | Return value of `prewarm_nbody_ephemeris_cache`. Returns `SimulationEngine.prewarm_nbody_ephemeris_cache(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.spaceagora|SpaceAGORA]] · `api` → `module_api` · call · `src/SpaceAGORA.jl`

**Downstream**

- `callees` → [[simulation.public_api_load_nbody_ephemeris_cache_bang|load_nbody_ephemeris_cache!]] · `callers` · call · `src/SpaceAGORA.jl:536-536`
- `callees` → [[spaceagora.spaceagora_load_nbody_ephemeris_cache_bang|load_nbody_ephemeris_cache!]] · `callers` · call · `src/SpaceAGORA.jl:536-536`
<!-- vulcan:connections:end -->

## Limitations
The cache is keyed on body set, start epoch, span and spacing; a later run with any different value silently misses the cache and rebuilds from SPICE. Registration is per Julia process, so multi-process campaigns must either call this on every worker or load a saved file. Prewarming a long mission at fine `dt_s` can allocate large arrays, and there is no size guard. The forwarding wrapper hides the keyword signature and provides no argument validation itself.

## Provenance
Mapped from `src/SpaceAGORA.jl` line 533.
