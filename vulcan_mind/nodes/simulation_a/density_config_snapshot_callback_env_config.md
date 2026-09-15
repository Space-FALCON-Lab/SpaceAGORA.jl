---
id: simulation_a.density_config_snapshot_callback_env_config
label: _snapshot_callback_env_config
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _snapshot_callback_env_config
  lines:
  - 177
  - 205
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
  description: Process environment variables holding the `SPACEAGORA_*` density, GRAM,
    thread-policy and cache knobs consulted at run setup.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: env_config
  type: CallbackEnvConfig
  units: n/a
  description: Immutable typed snapshot of every env-derived callback knob, stored
    on `SharedBuffers` and read by hot paths through `_callback_env_config(p)`.
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
# _snapshot_callback_env_config

## Purpose
`_snapshot_callback_env_config` resolves every environment-variable knob that callbacks and the right-hand-side atmosphere sampler would otherwise consult per invocation into a single typed `CallbackEnvConfig` record. It is evaluated once during `run_simulation` setup, and the resulting struct is stored on the run's shared buffers so hot paths read plain struct fields instead of re-parsing strings.

## Model & Assumptions
The snapshot is the mechanism that guarantees the density callback and the RHS-side atmosphere sampling observe identical settings for the whole run. Parsing at construction time instead of at call time also means a malformed value fails fast with a descriptive `ArgumentError` from `_parse_bool_env`, `_parse_float_env` or the integer parser, rather than surfacing mid-solve. Booleans accept `1/0`, `true/false`, `yes/no` and `on/off`; unrecognised values are rejected outright rather than silently defaulted.

## Design & Implementation
The function is a flat constructor call listing twenty-five fields in declaration order: the nested `GramTrackCacheConfig`, GRAM profiling and time-window flags, the per-step density freeze switch, the four vacuum-GRAM cache knobs, and matched parallel-mode, thread-threshold, allow-with-outer and assume-threadsafe triples for the density, control and thermal callbacks plus the density batch and GRAM isolated-pool settings. The companion accessor `_callback_env_config(p)` reads the snapshot off `p.shared_buffers` and falls back to a fresh snapshot when parameters were hand-constructed, which keeps unit tests and `withenv` probes working without a full engine setup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `environment` | ENV | n/a | yes | Process environment variables holding the `SPACEAGORA_*` density, GRAM, thread-policy and cache knobs consulted at run setup. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `env_config` | CallbackEnvConfig | n/a | — | Immutable typed snapshot of every env-derived callback knob, stored on `SharedBuffers` and read by hot paths through `_callback_env_config(p)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__callback_env_config|_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:217-217`
- [[simulation.config__gram_track_trajectory_supported|_gram_track_trajectory_supported]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:170-170`
- [[simulation.setup__initialize_runtime_env_config_bang|_initialize_runtime_env_config!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:915-915`

**Downstream**

- `callees` → [[core.runtime_types_callbackenvconfig|CallbackEnvConfig]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:178-178`
- `callees` → [[simulation.config__control_callback_allow_with_outer|_control_callback_allow_with_outer]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:199-199`
- `callees` → [[simulation.config__control_callback_parallel_mode|_control_callback_parallel_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:197-197`
- `callees` → [[simulation.config__control_callback_thread_threshold|_control_callback_thread_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:198-198`
- `callees` → [[simulation.config__density_batch_mode|_density_batch_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:192-192`
- `callees` → [[simulation.config__density_batch_threshold|_density_batch_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:193-193`
- `callees` → [[simulation.config__density_callback_allow_with_outer|_density_callback_allow_with_outer]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:190-190`
- `callees` → [[simulation.config__density_callback_parallel_mode|_density_callback_parallel_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:188-188`
- `callees` → [[simulation.config__density_callback_thread_threshold|_density_callback_thread_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:189-189`
- `callees` → [[simulation.config__density_freeze_per_step_enabled|_density_freeze_per_step_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:183-183`
- `callees` → [[simulation.config__gram_isolated_pool_max_workers|_gram_isolated_pool_max_workers]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:196-196`
- `callees` → [[simulation.config__gram_isolated_pool_mode|_gram_isolated_pool_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:194-194`
- `callees` → [[simulation.config__gram_isolated_pool_threshold|_gram_isolated_pool_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:195-195`
- `callees` → [[simulation.config__gram_track_cache_ignore_time_window|_gram_track_cache_ignore_time_window]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:181-181`
- `callees` → [[simulation.config__gram_track_cache_target_use_j2|_gram_track_cache_target_use_j2]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:182-182`
- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:191-191`
- `callees` → [[simulation.config__thermal_callback_allow_with_outer|_thermal_callback_allow_with_outer]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:203-203`
- `callees` → [[simulation.config__thermal_callback_parallel_mode|_thermal_callback_parallel_mode]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:201-201`
- `callees` → [[simulation.config__thermal_callback_thread_threshold|_thermal_callback_thread_threshold]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:202-202`
- `callees` → [[simulation.registry__gram_runtime_stats_enabled|_gram_runtime_stats_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:180-180`
- `callees` → [[simulation.vacuum_predicted_gram__vacuum_gram_cache_enabled|_vacuum_gram_cache_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:184-184`
- `callees` → [[simulation.vacuum_predicted_gram__vacuum_gram_cache_npoints|_vacuum_gram_cache_npoints]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:185-185`
- `callees` → [[simulation_a.gram_cache_config_gram_track_cache_config|_gram_track_cache_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:179-179`
<!-- vulcan:connections:end -->

## Limitations
Because the snapshot is taken once, changing an environment variable after a solve has started has no effect on that solve. The fallback path in `_callback_env_config` re-parses ENV on every call, so hand-constructed parameter objects lose the performance benefit the snapshot exists to provide. Field order in the constructor is positional and must stay synchronised with the `CallbackEnvConfig` definition in `runtime_types.jl`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl:177-205`.
