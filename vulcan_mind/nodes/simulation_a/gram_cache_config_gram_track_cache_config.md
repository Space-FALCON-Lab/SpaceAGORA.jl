---
id: simulation_a.gram_cache_config_gram_track_cache_config
label: _gram_track_cache_config
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/config.jl
  symbol: _gram_track_cache_config
  lines:
  - 69
  - 157
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: environment
  type: ENV
  units: n/a
  required: true
  description: Process environment holding the `SPACEAGORA_GRAM_TRACK_CACHE*` and
    legacy `SPACEAGORA_GRAM_SEGMENT_CACHE*` knobs for mode, horizons, tolerances and
    point counts.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: track_cache_config
  type: GramTrackCacheConfig
  units: n/a
  description: Typed configuration carrying the cache mode plus separate entry and
    orbit horizons, altitude tolerances, angular tolerances and point counts, and
    the transition band width.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# _gram_track_cache_config

## Purpose
`_gram_track_cache_config` resolves the GRAM track-cache settings from the environment into a typed `GramTrackCacheConfig`. It is called once during `_snapshot_callback_env_config`, so the cache behaviour observed by a solve is fixed at setup and cannot drift mid-run.

## Model & Assumptions
The configuration is regime-split: entry and orbit each get their own prediction horizon, altitude tolerance, angular tolerance and sample count, because a hypersonic entry pass needs a short, densely sampled segment while an orbital arc tolerates a long, sparse one. The defaults reflect that asymmetry — a one-second entry horizon with sixteen points against an eight-second orbital horizon with forty-eight points, and a five-hundred-metre entry altitude tolerance. The mode defaults to `:off` on the strength of a recorded benchmark in which a point-to-point entry case ran in roughly 0.45 s against roughly 42.8 s with the cache and Allen-Eggers targeting enabled.

## Design & Implementation
Legacy global knobs are honoured through a two-level `something` chain: a regime-specific variable wins, then the compatibility global, then the built-in default. Both the current `SPACEAGORA_GRAM_TRACK_CACHE` name and the older `SPACEAGORA_GRAM_SEGMENT_CACHE` name are accepted for the mode, which parses to `:off`, `:on` or `:auto` and throws `ArgumentError` on anything else. Numeric knobs go through `_parse_float_env_optional` and `_parse_int_env_optional`, which return `nothing` when unset so the fallback chain stays explicit, and every result is floored — horizons at one millisecond, point counts at two, tolerances at zero — so a hostile value cannot produce a degenerate cache. Angular tolerances are supplied in degrees and stored in radians.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `environment` | ENV | n/a | yes | Process environment holding the `SPACEAGORA_GRAM_TRACK_CACHE*` and legacy `SPACEAGORA_GRAM_SEGMENT_CACHE*` knobs for mode, horizons, tolerances and point counts. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `track_cache_config` | GramTrackCacheConfig | n/a | — | Typed configuration carrying the cache mode plus separate entry and orbit horizons, altitude tolerances, angular tolerances and point counts, and the transition band width. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:179-179`

**Downstream**

- `callees` → [[core.runtime_types_gramtrackcacheconfig|GramTrackCacheConfig]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:145-145`
- `callees` → [[simulation.config__gram_track_cache_mode|_gram_track_cache_mode]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:70-70`
- `callees` → [[simulation.config__parse_float_env|_parse_float_env]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:143-143`
- `callees` → [[simulation.config__parse_float_env_optional|_parse_float_env_optional]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:73-73`
- `callees` → [[simulation.config__parse_int_env_optional|_parse_int_env_optional]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:76-76`
<!-- vulcan:connections:end -->

## Limitations
Because the values are snapshotted, changing an environment variable after setup has no effect on the running solve. The compatibility chain means a legacy global silently overrides the built-in default for both regimes at once, which can surprise a caller who set only the entry-specific variable expecting the orbit default to remain. The function validates ranges but not physical consistency between the entry and orbit regimes.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/config.jl:68-157`.
