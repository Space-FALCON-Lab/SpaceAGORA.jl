---
id: core.runtime_types_callbackenvconfig
label: CallbackEnvConfig
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: CallbackEnvConfig
  lines:
  - 605
  - 605
inputs:
- id: gram_track_cache
  type: GramTrackCacheConfig
  units: n/a
  required: true
  description: Field `gram_track_cache`.
- id: gram_runtime_stats_enabled
  type: Bool
  units: n/a
  required: true
  description: Field `gram_runtime_stats_enabled`.
- id: gram_track_cache_ignore_time_window
  type: Bool
  units: n/a
  required: true
  description: Field `gram_track_cache_ignore_time_window`.
- id: gram_track_cache_target_use_j2
  type: Bool
  units: n/a
  required: true
  description: Field `gram_track_cache_target_use_j2`.
- id: density_freeze_per_step
  type: Bool
  units: n/a
  required: true
  description: Field `density_freeze_per_step`.
- id: vacuum_gram_cache_enabled
  type: Bool
  units: n/a
  required: true
  description: Field `vacuum_gram_cache_enabled`.
- id: vacuum_gram_cache_npoints
  type: Int
  units: n/a
  required: true
  description: Field `vacuum_gram_cache_npoints`.
- id: vacuum_gram_cache_horizon_s
  type: Float64
  units: n/a
  required: true
  description: Field `vacuum_gram_cache_horizon_s`.
- id: vacuum_gram_cache_deviation_m
  type: Float64
  units: n/a
  required: true
  description: Field `vacuum_gram_cache_deviation_m`.
- id: density_parallel_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `density_parallel_mode`.
- id: density_thread_threshold
  type: Int
  units: n/a
  required: true
  description: Field `density_thread_threshold`.
- id: density_allow_with_outer
  type: Bool
  units: n/a
  required: true
  description: Field `density_allow_with_outer`.
- id: density_assume_threadsafe
  type: Bool
  units: n/a
  required: true
  description: Field `density_assume_threadsafe`.
- id: density_batch_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `density_batch_mode`.
- id: density_batch_threshold
  type: Int
  units: n/a
  required: true
  description: Field `density_batch_threshold`.
- id: gram_isolated_pool_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `gram_isolated_pool_mode`.
- id: gram_isolated_pool_threshold
  type: Int
  units: n/a
  required: true
  description: Field `gram_isolated_pool_threshold`.
- id: gram_isolated_pool_max_workers
  type: Int
  units: n/a
  required: true
  description: Field `gram_isolated_pool_max_workers`.
- id: control_parallel_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `control_parallel_mode`.
- id: control_thread_threshold
  type: Int
  units: n/a
  required: true
  description: Field `control_thread_threshold`.
- id: control_allow_with_outer
  type: Bool
  units: n/a
  required: true
  description: Field `control_allow_with_outer`.
- id: control_assume_threadsafe
  type: Bool
  units: n/a
  required: true
  description: Field `control_assume_threadsafe`.
- id: thermal_parallel_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `thermal_parallel_mode`.
- id: thermal_thread_threshold
  type: Int
  units: n/a
  required: true
  description: Field `thermal_thread_threshold`.
- id: thermal_allow_with_outer
  type: Bool
  units: n/a
  required: true
  description: Field `thermal_allow_with_outer`.
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
  type: CallbackEnvConfig
  units: n/a
  description: Constructed `CallbackEnvConfig`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# CallbackEnvConfig

## Purpose
Immutable run-scoped snapshot of every `SPACEAGORA_*` environment knob consulted by the density, control, and thermal callbacks and by RHS-side atmosphere sampling, so hot paths read struct fields instead of parsing `ENV`.

## Design & Implementation
Plain `struct CallbackEnvConfig` with a nested `gram_track_cache::GramTrackCacheConfig` and about thirty scalar fields: booleans such as `gram_runtime_stats_enabled`, `density_freeze_per_step` (reuse one density sample per accepted step across RK stages), `vacuum_gram_cache_enabled`, and `density_assume_threadsafe`; integers like `vacuum_gram_cache_npoints`, `density_thread_threshold`, `gram_isolated_pool_max_workers`; floats `vacuum_gram_cache_horizon_s` and `vacuum_gram_cache_deviation_m` (m); and `Symbol` modes `density_parallel_mode`, `density_batch_mode`, `gram_isolated_pool_mode`, `control_parallel_mode`, `thermal_parallel_mode`. Built once by `SimulationCallbacks` at `run_simulation` setup and stored in `SharedBuffers.callback_env_config[]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `gram_track_cache` | GramTrackCacheConfig | n/a | yes | Field `gram_track_cache`. |
| in | `gram_runtime_stats_enabled` | Bool | n/a | yes | Field `gram_runtime_stats_enabled`. |
| in | `gram_track_cache_ignore_time_window` | Bool | n/a | yes | Field `gram_track_cache_ignore_time_window`. |
| in | `gram_track_cache_target_use_j2` | Bool | n/a | yes | Field `gram_track_cache_target_use_j2`. |
| in | `density_freeze_per_step` | Bool | n/a | yes | Field `density_freeze_per_step`. |
| in | `vacuum_gram_cache_enabled` | Bool | n/a | yes | Field `vacuum_gram_cache_enabled`. |
| in | `vacuum_gram_cache_npoints` | Int | n/a | yes | Field `vacuum_gram_cache_npoints`. |
| in | `vacuum_gram_cache_horizon_s` | Float64 | n/a | yes | Field `vacuum_gram_cache_horizon_s`. |
| in | `vacuum_gram_cache_deviation_m` | Float64 | n/a | yes | Field `vacuum_gram_cache_deviation_m`. |
| in | `density_parallel_mode` | Symbol | n/a | yes | Field `density_parallel_mode`. |
| in | `density_thread_threshold` | Int | n/a | yes | Field `density_thread_threshold`. |
| in | `density_allow_with_outer` | Bool | n/a | yes | Field `density_allow_with_outer`. |
| in | `density_assume_threadsafe` | Bool | n/a | yes | Field `density_assume_threadsafe`. |
| in | `density_batch_mode` | Symbol | n/a | yes | Field `density_batch_mode`. |
| in | `density_batch_threshold` | Int | n/a | yes | Field `density_batch_threshold`. |
| in | `gram_isolated_pool_mode` | Symbol | n/a | yes | Field `gram_isolated_pool_mode`. |
| in | `gram_isolated_pool_threshold` | Int | n/a | yes | Field `gram_isolated_pool_threshold`. |
| in | `gram_isolated_pool_max_workers` | Int | n/a | yes | Field `gram_isolated_pool_max_workers`. |
| in | `control_parallel_mode` | Symbol | n/a | yes | Field `control_parallel_mode`. |
| in | `control_thread_threshold` | Int | n/a | yes | Field `control_thread_threshold`. |
| in | `control_allow_with_outer` | Bool | n/a | yes | Field `control_allow_with_outer`. |
| in | `control_assume_threadsafe` | Bool | n/a | yes | Field `control_assume_threadsafe`. |
| in | `thermal_parallel_mode` | Symbol | n/a | yes | Field `thermal_parallel_mode`. |
| in | `thermal_thread_threshold` | Int | n/a | yes | Field `thermal_thread_threshold`. |
| in | `thermal_allow_with_outer` | Bool | n/a | yes | Field `thermal_allow_with_outer`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CallbackEnvConfig | n/a | — | Constructed `CallbackEnvConfig`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:178-178`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:178-178`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Values are captured at run start, so changing `ENV` during a run has no effect; accessor code falls back to live parsing only when the ref is still `nothing`. There are no field-level invariants (thresholds may be negative, modes unknown), validation being the builder's responsibility. Adding a knob requires editing this struct, its builder, and every accessor.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 605.
