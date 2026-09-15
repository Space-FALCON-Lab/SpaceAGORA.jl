---
id: simulation.runtime__density_state_from_kinematics_bang
label: _density_state_from_kinematics!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: _density_state_from_kinematics!
  lines:
  - 78
  - 78
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: vel_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel_ii`.
- id: current_mass_kg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `current_mass_kg`.
- id: alt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt`.
- id: lat
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat`.
- id: lon
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: density_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `density_model`.
- id: cache_cfg
  type: Any
  units: n/a
  required: true
  description: Positional argument `cache_cfg`.
- id: stats_enabled
  type: Bool
  units: n/a
  required: true
  description: Positional argument `stats_enabled`.
- id: target_include_j2
  type: Bool
  units: n/a
  required: true
  description: Positional argument `target_include_j2`.
- id: caches
  type: Vector{Union{Nothing, GramTrackCache}}
  units: n/a
  required: true
  description: Positional argument `caches`.
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_density_state_from_kinematics!`; mutates `p` in place.
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

# _density_state_from_kinematics!

## Purpose
The core density decision for one satellite: answer from the vacuum-predicted cache, the GRAM track cache, or a direct model call, updating runtime statistics as it goes.

## Design & Implementation
Three tiers. If the run-scoped `vacuum_gram_cache_enabled` flag is set and the satellite's `in_atmosphere` flag is true, it queries the log-density spline cache and returns. Otherwise, if the track cache applies to this density model, it either probes the cache for a segment bracketing `t` — with a stats path that additionally measures altitude, latitude and longitude interpolation error and classifies hits and misses by time window or state tolerance — or uses the lighter `_gram_track_cache_ready` check. A usable segment is evaluated by `_gram_track_cache_eval`; otherwise `_gram_track_cache_refresh!` rebuilds the cache out to `_density_segment_end_t` and returns the fresh value. With no cache the model's `getDensity` is called directly. Stats updates go through `_gram_runtime_stats_update!` closures.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `vel_ii` | SVector{3, Float64} | n/a | yes | Positional argument `vel_ii`. |
| in | `current_mass_kg` | Float64 | n/a | yes | Positional argument `current_mass_kg`. |
| in | `alt` | Float64 | n/a | yes | Positional argument `alt`. |
| in | `lat` | Float64 | n/a | yes | Positional argument `lat`. |
| in | `lon` | Float64 | n/a | yes | Positional argument `lon`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `density_model` | Any | n/a | yes | Positional argument `density_model`. |
| in | `cache_cfg` | Any | n/a | yes | Positional argument `cache_cfg`. |
| in | `stats_enabled` | Bool | n/a | yes | Positional argument `stats_enabled`. |
| in | `target_include_j2` | Bool | n/a | yes | Positional argument `target_include_j2`. |
| in | `caches` | Vector{Union{Nothing, GramTrackCache}} | n/a | yes | Positional argument `caches`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_density_state_from_kinematics!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:227-227`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:227-227`

**Downstream**

- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:201-201`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:94-94`
- `callees` → [[simulation.interpolation__angdiff_rad|_angdiff_rad]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:140-140`
- `callees` → [[simulation.interpolation__gram_track_cache_enabled|_gram_track_cache_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:118-118`
- `callees` → [[simulation.interpolation__gram_track_cache_eval|_gram_track_cache_eval]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:171-171`
- `callees` → [[simulation.interpolation__gram_track_cache_profile|_gram_track_cache_profile]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:124-124`
- `callees` → [[simulation.interpolation__gram_track_cache_segment|_gram_track_cache_segment]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:127-127`
- `callees` → [[simulation.interpolation__lerp|_lerp]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:136-136`
- `callees` → [[simulation.interpolation__lerp_angle_rad|_lerp_angle_rad]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:137-137`
- `callees` → [[simulation.model_selection__gram_density_cache_for_sat_bang|_gram_density_cache_for_sat!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:125-125`
- `callees` → [[simulation.registry__gram_runtime_stats_update_bang|_gram_runtime_stats_update!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:113-113`
- `callees` → [[simulation.runtime__density_segment_end_t|_density_segment_end_t]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:173-173`
- `callees` → [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:103-103`
- `callees` → [[simulation.vacuum_predicted_gram__vacuum_gram_cache_for_sat_bang|_vacuum_gram_cache_for_sat!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:102-102`
- `callees` → [[simulation_a.interpolation_gram_track_cache_ready|_gram_track_cache_ready]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:167-167`
- `callees` → [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:174-174`
<!-- vulcan:connections:end -->

## Limitations
The vacuum cache is consulted only inside the atmosphere, so the `in_atmosphere` flag must already be correct when this runs; the stats-enabled probe duplicates the tolerance logic of `_gram_track_cache_ready`, so the two can drift apart. The function mutates the per-satellite cache entries in `caches` and the runtime stats as side effects.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 78.
