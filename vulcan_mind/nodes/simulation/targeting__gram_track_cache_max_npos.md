---
id: simulation.targeting__gram_track_cache_max_npos
label: _gram_track_cache_max_npos
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_track_cache_max_npos
  lines:
  - 3
  - 3
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
  description: Return value of `_gram_track_cache_max_npos`.
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

# _gram_track_cache_max_npos

## Purpose
Reads the upper bound on the number of positions per GRAM track-cache request from the environment variable `SPACEAGORA_GRAM_TRACK_CACHE_MAX_NPOS`, defaulting to 512.

## Design & Implementation
Fetches and `strip`s the variable with `get(ENV, ..., "512")`, then `parse(Int, raw)` inside a `try`. A parse failure is converted into an `ArgumentError` whose message includes the offending string. The parsed value is floored with `max(2, parsed)` so at least two positions are always allowed (a track needs a start and an end). The function is `@inline` and returns `Int`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_gram_track_cache_max_npos`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:174-174`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The environment is consulted on every call, so this should not be invoked in a hot loop; callers cache the result in `CallbackEnvConfig`. Negative or zero values are silently promoted to 2 instead of being rejected. There is no upper bound, so an enormous value can cause the track cache to request far more GRAM evaluations than memory or time allow. Non-ASCII whitespace is handled by `strip`, but embedded plus signs or underscores are rejected by `parse`.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 3.
