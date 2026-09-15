---
id: simulation.registry__gram_runtime_stats_enabled
label: _gram_runtime_stats_enabled
kind: function
source:
  file: src/simulation/callbacks/registry.jl
  symbol: _gram_runtime_stats_enabled
  lines:
  - 65
  - 65
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
  description: Return value of `_gram_runtime_stats_enabled`.
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

# _gram_runtime_stats_enabled

## Purpose
Gate that decides whether the GRAM atmosphere runtime profiler collects statistics during a simulation, read on each query from the `SPACEAGORA_GRAM_PROFILE` environment variable.

## Design & Implementation
An `@inline` zero-argument function returning `Bool`. It delegates to `_parse_bool_env("SPACEAGORA_GRAM_PROFILE", false)`, so an unset or unparseable value yields `false` and profiling stays off. Callers in the density callbacks wrap `_gram_runtime_stats_update!` calls in this check so the `ReentrantLock` in `_gram_runtime_stats_lock` is never taken on the hot path unless profiling is explicitly requested.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gram_runtime_stats_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:180-180`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:77-77`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/registry.jl:66-66`
<!-- vulcan:connections:end -->

## Limitations
The environment variable is re-read on every call, so `ENV` lookups occur at each density evaluation rather than being cached at simulation start; flipping the variable mid-run changes behaviour immediately. The variable name is hard-coded and the accepted truthy spellings are whatever `_parse_bool_env` recognises.

## Provenance
Mapped from `src/simulation/callbacks/registry.jl` line 65.
