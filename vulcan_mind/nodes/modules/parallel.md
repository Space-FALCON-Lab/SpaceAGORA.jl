---
id: module.parallel
label: ParallelProfiles
kind: module
source:
  file: src/parallel/routing/parallel_profiles.jl
  symbol: ParallelProfiles
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: Parallel profile configuration, environment mapping, outer-route state
    and metrics, route selection, and process-pool APIs.
tags:
- module
charts:
- master
origin: agent
---

# ParallelProfiles

## Purpose
`ParallelProfiles` is the configuration and routing boundary for process- and thread-aware campaign execution. Its aggregator includes profile definitions, environment mapping, outer-route state, route selection, and route metrics, then exports their public records and functions. Simulation campaigns use this namespace to choose an execution route without embedding environment-variable parsing or process-pool lifecycle rules in the campaign algorithms.

## Theory & Math
Route selection scores candidate execution routes from observed features and tuning parameters. The resulting route is discrete, while metrics such as throughput and failure counts are accumulated over campaign calls. Profile configuration is loaded from TOML or environment pairs; the numerical simulation itself remains unchanged by the route choice.

## Model & Assumptions
The profile name and environment mapping must be internally consistent. Outer-route state is mutable process-local state and must be reset or loaded deliberately between campaigns. A process pool is assumed to be safe for the requested worker count and compatible with the Julia serialization boundary used by the campaign function.

## Design & Implementation
`parallel_profiles.jl` includes five implementation files. `ParallelProfile` and `profile_config` define and resolve profile data. `default_outer_route` and `select_outer_route!` choose a route from candidate state. `ProcessPool` and `ensure_process_workers!` manage worker availability. The exported `with_parallel_profile` scopes environment changes, while route statistics and persistence functions expose observability and restart support.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | Parallel profile configuration, environment mapping, outer-route state and metrics, route selection, and process-pool APIs. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[grp.src_parallel_policy|parallel/policy/]] · `members_in` · call · `src/parallel/policy/adaptive_decision.jl`
- `api` → [[grp.src_parallel_process|parallel/process/]] · `members_in` · call · `src/parallel/process/worker_pool.jl`
- `api` → [[grp.src_parallel_routing|parallel/routing/]] · `members_in` · call · `src/parallel/routing/env_mapping.jl`
- `api` → [[module.spaceagora|SpaceAGORA]] · `parallel` · call · `src/SpaceAGORA.jl:6-7`
- `api` → [[parallel.adaptive_decision_use_threads_policy|use_threads_policy]] · `module_api` · call · `src/parallel/policy/adaptive_decision.jl`
- `api` → [[parallel.context__active_policy_scope_id|_active_policy_scope_id]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__create_persistent_foreach_pool|_create_persistent_foreach_pool]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__create_persistent_foreach_worker_pool|_create_persistent_foreach_worker_pool]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__create_spin_barrier_pool|_create_spin_barrier_pool]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__persistent_foreach_worker_loop|_persistent_foreach_worker_loop]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__persistent_foreach_worker_loop_w|_persistent_foreach_worker_loop_w]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__shutdown_persistent_foreach_pool_bang|_shutdown_persistent_foreach_pool!]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__shutdown_spin_barrier_pool_bang|_shutdown_spin_barrier_pool!]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__spin_barrier_dispatch_bang|_spin_barrier_dispatch!]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__spin_barrier_pool_for|_spin_barrier_pool_for]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.context__spin_barrier_worker_loop_w|_spin_barrier_worker_loop_w]] · `module_api` · call · `src/parallel/policy/context.jl`
- `api` → [[parallel.env_config__default_thread_pool_size|_default_thread_pool_size]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_config__persistent_hint_default_path|_persistent_hint_default_path]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_config__persistent_hint_path|_persistent_hint_path]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_config__safe_token|_safe_token]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_config__telemetry_bucket|_telemetry_bucket]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_config_persistent_hints_persist_enabled|persistent_hints_persist_enabled]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_config_persistent_hints_state_reset_requested|persistent_hints_state_reset_requested]] · `module_api` · call · `src/parallel/policy/env_config.jl`
- `api` → [[parallel.env_mapping__machine_parallel_class|_machine_parallel_class]] · `module_api` · call · `src/parallel/routing/env_mapping.jl`
- `api` → [[parallel.observation_tracking_record_route_discard_bang|record_route_discard!]] · `module_api` · call · `src/parallel/policy/observation_tracking.jl`
- `api` → [[parallel.outer_route_selection__best_candidate|_best_candidate]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__candidate_confidence_width|_candidate_confidence_width]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__compat_outer_route_signature|_compat_outer_route_signature]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__feature_heavy_for_process|_feature_heavy_for_process]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_count_bucket|_route_count_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_density_bucket|_route_density_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_effector_cost_bucket|_route_effector_cost_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_elapsed_stats|_route_elapsed_stats]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_harmonics_bucket|_route_harmonics_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_interval_bucket|_route_interval_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_link_bucket|_route_link_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_max_link_bucket|_route_max_link_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_mission_bucket|_route_mission_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_ranked_candidates|_route_ranked_candidates]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_sat_bucket|_route_sat_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection__route_solver_bucket|_route_solver_bucket]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_selection_outer_route_stats_snapshot|outer_route_stats_snapshot]] · `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- `api` → [[parallel.outer_route_state__route_payload_stats|_route_payload_stats]] · `module_api` · call · `src/parallel/routing/outer_route_state.jl`
- `api` → [[parallel.outer_route_state__route_stats_payload|_route_stats_payload]] · `module_api` · call · `src/parallel/routing/outer_route_state.jl`
- `api` → [[parallel.outer_route_state_load_outer_route_state_bang|load_outer_route_state!]] · `module_api` · call · `src/parallel/routing/outer_route_state.jl`
- `api` → [[parallel.outer_route_state_save_outer_route_state|save_outer_route_state]] · `module_api` · call · `src/parallel/routing/outer_route_state.jl`
- `api` → [[parallel.parallel_profile|ParallelProfile]] · `module_api` · call · `src/parallel/routing/profile_definitions.jl`
- `api` → [[parallel.persistent_hints__hint_bucket|_hint_bucket]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.persistent_hints__hint_payload_stats|_hint_payload_stats]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.persistent_hints__hint_signature_value|_hint_signature_value]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.persistent_hints__hint_stats_payload|_hint_stats_payload]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.persistent_hints__load_persistent_hint_state_locked_bang|_load_persistent_hint_state_locked!]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.persistent_hints__save_persistent_hint_state_locked_bang|_save_persistent_hint_state_locked!]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.persistent_hints_reset_persistent_hint_state_bang|reset_persistent_hint_state!]] · `module_api` · call · `src/parallel/policy/persistent_hints.jl`
- `api` → [[parallel.policy_telemetry_reset_policy_telemetry_bang|reset_policy_telemetry!]] · `module_api` · call · `src/parallel/policy/policy_telemetry.jl`
- `api` → [[parallel.process_pool|ProcessPool]] · `module_api` · call · `src/parallel/process/worker_pool.jl`
- `api` → [[parallel.profile_definitions__normalize_profile_token|_normalize_profile_token]] · `module_api` · call · `src/parallel/routing/profile_definitions.jl`
- `api` → [[parallel.thread_execution__persistent_pool_for|_persistent_pool_for]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution__persistent_pool_key|_persistent_pool_key]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution__persistent_worker_pool_for|_persistent_worker_pool_for]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution__threaded_foreach_persistent_bang|_threaded_foreach_persistent!]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_thread_worker_count|thread_worker_count]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_threaded_collect_bang|threaded_collect!]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_threaded_collect_persistent_bang|threaded_collect_persistent!]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_threaded_foreach_worker_spin|threaded_foreach_worker_spin]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.thread_execution_threaded_reduce|threaded_reduce]] · `module_api` · call · `src/parallel/policy/thread_execution.jl`
- `api` → [[parallel.types__hintlayerstatsaccumulator|_HintLayerStatsAccumulator]] · `module_api` · call · `src/parallel/policy/types.jl`
- `api` → [[parallel.types__persistentforeachpool|_PersistentForeachPool]] · `module_api` · call · `src/parallel/policy/types.jl`
- `api` → [[parallel.types__persistenthintstate|_PersistentHintState]] · `module_api` · call · `src/parallel/policy/types.jl`
- `api` → [[parallel.types__spinbarrierpool|_SpinBarrierPool]] · `module_api` · call · `src/parallel/policy/types.jl`
- `api` → [[parallel.types_adaptivechoicestats|AdaptiveChoiceStats]] · `module_api` · call · `src/parallel/policy/types.jl`
- `api` → [[parallel.types_adaptivecontrollerstate|AdaptiveControllerState]] · `module_api` · call · `src/parallel/policy/types.jl`
- `api` → [[parallel.worker_pool__furnish_default_spice_kernels_bang|_furnish_default_spice_kernels!]] · `module_api` · call · `src/parallel/process/worker_pool.jl`
- `api` → [[parallel.worker_pool__warm_gram_wrapper_bang|_warm_gram_wrapper!]] · `module_api` · call · `src/parallel/process/worker_pool.jl`
- `api` → [[parallel.worker_pool_shutdown_process_pool_bang|shutdown_process_pool!]] · `module_api` · call · `src/parallel/process/worker_pool.jl`
- `api` → [[parcore.parallel_policy_parallelpolicy|ParallelPolicy]] · `module_api` · call · `src/parallel/policy/parallel_policy.jl`
- `api` → [[parcore.parallel_process_parallelprocess|ParallelProcess]] · `module_api` · call · `src/parallel/process/parallel_process.jl`
- `api` → [[parcore.parallel_profiles_parallelprofiles|ParallelProfiles]] · `module_api` · call · `src/parallel/routing/parallel_profiles.jl`
- `api` → [[parcore.types_policytelemetry|PolicyTelemetry]] · `module_api` · call · `src/parallel/policy/types.jl`
<!-- vulcan:connections:end -->

## Limitations
Process startup cost, worker availability, environment inheritance, and serialized closure size can dominate campaign runtime. Outer-route feedback is only as meaningful as the features recorded by the caller. A route selected from stale persisted state can be inappropriate for a new machine or workload, so callers should load state with an explicit compatibility policy.

## Provenance
Mapped from `src/parallel/routing/parallel_profiles.jl` and its included routing and process implementation files.
