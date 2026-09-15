---
id: simulation.config__parse_bool_env
label: _parse_bool_env
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _parse_bool_env
  lines:
  - 1
  - 1
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Bool
  units: n/a
  required: true
  description: Positional argument `default`.
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
  description: Return value of `_parse_bool_env`.
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

# _parse_bool_env

## Purpose
Single boolean parser for every `SPACEAGORA_*` environment switch consulted by the density/control/thermal callback layer, so that all knobs accept the same spelling set and fail the same way.

## Design & Implementation
Reads `ENV[name]`, substituting the string `"1"` or `"0"` derived from `default` when the key is absent, then lowercases and strips it. Accepts `1/true/yes/on` as true and `0/false/no/off` as false. Anything else raises `ArgumentError` naming the variable and the offending raw value. Marked `@inline` because callers invoke it inside `_snapshot_callback_env_config` and, on fallback paths, per callback invocation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Bool | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_parse_bool_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:39-39`
- [[simulation.config__control_callback_allow_with_outer|_control_callback_allow_with_outer]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:107-107`
- [[simulation.config__density_callback_allow_with_outer|_density_callback_allow_with_outer]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:55-55`
- [[simulation.config__density_freeze_per_step_enabled|_density_freeze_per_step_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:21-21`
- [[simulation.config__gram_track_cache_ignore_time_window|_gram_track_cache_ignore_time_window]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:12-12`
- [[simulation.config__gram_track_cache_target_use_j2|_gram_track_cache_target_use_j2]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:16-16`
- [[simulation.config__thermal_callback_allow_with_outer|_thermal_callback_allow_with_outer]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:126-126`
- [[simulation.registry__gram_runtime_stats_enabled|_gram_runtime_stats_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/registry.jl:66-66`
- [[simulation.targeting__gram_track_cache_periapsis_split_enabled|_gram_track_cache_periapsis_split_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:1-1`
- [[simulation.vacuum_predicted_gram__vacuum_gram_cache_enabled|_vacuum_gram_cache_enabled]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:19-19`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:191-191`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Whitespace-only or mixed-case input is normalised, but any other token throws rather than falling back to the default, so a typo in a shell export aborts the run at snapshot time instead of degrading silently. Reading `ENV` is not synchronised, so a concurrent `withenv` in another task can be observed mid-change.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 1.
