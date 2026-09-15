---
id: simulation.public_api_prewarm_nbody_ephemeris_cache
label: prewarm_nbody_ephemeris_cache
kind: function
source:
  file: src/simulation/engine/public_api.jl
  symbol: prewarm_nbody_ephemeris_cache
  lines:
  - 35
  - 35
inputs:
- id: config
  type: SimulationEngineConfig
  units: n/a
  required: true
  description: Positional argument `config`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `prewarm_nbody_ephemeris_cache`. Returns `_with_engine_env_overrides(config,
    () -> prewarm_nbody_ephemeris_cache(args; kwa`.
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

# prewarm_nbody_ephemeris_cache

## Purpose
Public entry point for populating the n-body ephemeris cache ahead of a run, so that gravitational third-body lookups during integration hit memory instead of the underlying ephemeris kernels.

## Design & Implementation
Two methods are provided. The first takes a `SimulationEngineConfig` and forwards through `_with_engine_env_overrides(config, () -> prewarm_nbody_ephemeris_cache(args; kwargs...))`, so engine-level environment settings are in force for the duration of the prewarm and restored afterwards. The second is the working method, with keyword arguments `dt_s::Union{Nothing, Real}` for the sampling step in seconds, `mission_end_s::Union{Nothing, Real}` for the horizon in seconds, and `save_path::Union{Nothing, AbstractString}` for an optional cache file; all three default to `nothing` and are passed straight through to the internal `_prewarm_nbody_ephemeris_cache`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | SimulationEngineConfig | n/a | yes | Positional argument `config`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `prewarm_nbody_ephemeris_cache`. Returns `_with_engine_env_overrides(config, () -> prewarm_nbody_ephemeris_cache(args; kwa`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/public_api.jl`

**Downstream**

- `callees` → [[simulation.from_env__with_engine_env_overrides|_with_engine_env_overrides]] · `callers` · call · `src/simulation/engine/public_api.jl:36-36`
- `callees` → [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callers` · call · `src/simulation/engine/public_api.jl:45-45`
<!-- vulcan:connections:end -->

## Limitations
The public wrapper performs no validation: a negative or zero `dt_s`, a `mission_end_s` earlier than the start epoch, or an unwritable `save_path` are only detected inside the internal implementation. The untyped `args` positional accepts any object, so a wrong configuration type is not caught at this boundary. The prewarm cost scales with `mission_end_s / dt_s` and can dominate short runs.

## Provenance
Mapped from `src/simulation/engine/public_api.jl` line 35.
