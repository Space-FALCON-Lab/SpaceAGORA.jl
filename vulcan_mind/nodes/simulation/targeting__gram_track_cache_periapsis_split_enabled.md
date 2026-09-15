---
id: simulation.targeting__gram_track_cache_periapsis_split_enabled
label: _gram_track_cache_periapsis_split_enabled
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_track_cache_periapsis_split_enabled
  lines:
  - 1
  - 1
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
  type: Any
  units: n/a
  description: Return value of `_gram_track_cache_periapsis_split_enabled`. Returns
    `_parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT", true)`.
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

# _gram_track_cache_periapsis_split_enabled

## Purpose
Boolean switch, read from `SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT`, that controls whether a GRAM track-cache request spanning periapsis is split into two segments so the dense low-altitude region receives finer sampling.

## Design & Implementation
An `@inline` one-liner delegating to `_parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT", true)`, which accepts the usual truthy/falsy spellings and returns the default `true` when the variable is unset or blank. No caching is performed at this level.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_track_cache_periapsis_split_enabled`. Returns `_parse_bool_env("SPACEAGORA_GRAM_TRACK_CACHE_PERIAPSIS_SPLIT", true)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:91-91`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:1-1`
<!-- vulcan:connections:end -->

## Limitations
Each call performs an `ENV` lookup and string parse; the callback layer is expected to snapshot the value into its environment config rather than call this per step. The default of `true` changes cache behaviour for all missions including circular orbits where the split provides no benefit and merely doubles the number of requests. Unrecognised values are handled by `_parse_bool_env`, whose fallback policy this function inherits without documenting.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 1.
