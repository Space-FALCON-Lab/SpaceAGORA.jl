---
id: simulation.config__gram_track_cache_ignore_time_window
label: _gram_track_cache_ignore_time_window
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_track_cache_ignore_time_window
  lines:
  - 11
  - 11
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
  description: Return value of `_gram_track_cache_ignore_time_window`.
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

# _gram_track_cache_ignore_time_window

## Purpose
Reports whether the GRAM track cache is allowed to reuse a cached trajectory sample regardless of how far the sample's timestamp is from the query time.

## Design & Implementation
Thin `@inline` wrapper over `_parse_bool_env` for `SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW`, defaulting to `true`. The value is captured once into the `CallbackEnvConfig` snapshot built by `_snapshot_callback_env_config`, so the hot path reads a struct field rather than touching `ENV`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gram_track_cache_ignore_time_window`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.interpolation__gram_track_cache_segment|_gram_track_cache_segment]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:52-52`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:181-181`
- [[simulation_a.interpolation_gram_track_cache_ready|_gram_track_cache_ready]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:100-100`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:12-12`
<!-- vulcan:connections:end -->

## Limitations
Defaulting to `true` means the time-window guard is off unless deliberately enabled, which can return atmospheric samples generated for a noticeably different epoch. The function performs no consistency check against the cache horizon settings, so contradictory combinations of cache knobs are accepted.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 11.
