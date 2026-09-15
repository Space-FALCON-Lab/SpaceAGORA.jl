---
id: parallel.env_config_parse_thread_threshold_env
label: parse_thread_threshold_env
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: parse_thread_threshold_env
  lines:
  - 23
  - 23
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
- id: default
  type: Int
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
  type: Int
  units: n/a
  description: Return value of `parse_thread_threshold_env`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# parse_thread_threshold_env

## Purpose
Parses an integer tuning knob (thread counts, minimum satellite counts, chunk sizes) from the environment and clamps it to at least 1 so downstream divisions and loop bounds never see zero.

## Design & Implementation
Takes `name::String` and `default::Int`; the default is stringified for `get(ENV, ...)` then stripped. `parse(Int, raw)` runs inside a `try` block whose `catch` rethrows as `ArgumentError("<name> must be an integer, got '<raw>'")`. The final `max(1, threshold)` silently lifts zero and negative values to 1. Used by `inner_dynamic_chunk_size`, `persistent_hints_min_samples`, `auto_thread_min_budget`, and `adaptive_window_size`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `default` | Int | n/a | yes | Positional argument `default`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `parse_thread_threshold_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_max_threads|_multibody_max_threads]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:17-17`
- [[dynamics.aerodynamic_wrench_models__multibody_thread_threshold|_multibody_thread_threshold]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:13-13`
- [[parallel.env_config__snapshot_auto_min_budget|_snapshot_auto_min_budget]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:209-209`
- [[parallel.env_config_auto_thread_min_budget|auto_thread_min_budget]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:137-137`
- [[parallel.env_config_inner_dynamic_chunk_size|inner_dynamic_chunk_size]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:64-64`
- [[parallel.env_config_persistent_hints_min_samples|persistent_hints_min_samples]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:89-89`
- [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:209-209`
- [[simulation.config__control_callback_thread_threshold|_control_callback_thread_threshold]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:103-103`
- [[simulation.config__density_batch_threshold|_density_batch_threshold]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:63-63`
- [[simulation.config__density_callback_thread_threshold|_density_callback_thread_threshold]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:51-51`
- [[simulation.config__gram_isolated_pool_max_workers|_gram_isolated_pool_max_workers]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:85-85`
- [[simulation.config__gram_isolated_pool_threshold|_gram_isolated_pool_threshold]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:81-81`
- [[simulation.config__thermal_callback_thread_threshold|_thermal_callback_thread_threshold]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:119-119`
- [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callees` → `callers` · call · `src/simulation/engine/adapters/from_env.jl:126-126`
- [[simulation.setup__effector_cost_min_samples|_effector_cost_min_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:422-422`
- [[simulation.setup__effector_long_orbit_threshold|_effector_long_orbit_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:474-474`
- [[simulation.setup__effector_max_threads|_effector_max_threads]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:402-402`
- [[simulation.setup__effector_thread_threshold|_effector_thread_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:398-398`
- [[simulation.setup__nbody_ephemeris_cache_max_samples|_nbody_ephemeris_cache_max_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:213-213`
- [[simulation.setup__planet_frame_cache_max_samples|_planet_frame_cache_max_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:225-225`
- [[simulation.setup__rhs_batch_thread_threshold|_rhs_batch_thread_threshold]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:372-372`
- [[simulation.setup__rhs_effector_cost_min_samples|_rhs_effector_cost_min_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:458-458`
- [[simulation.setup__rhs_flat_min_effectors|_rhs_flat_min_effectors]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:781-781`
- [[simulation.setup__rhs_flat_min_sats|_rhs_flat_min_sats]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:777-777`
- [[simulation.setup__rhs_flat_min_thread_budget|_rhs_flat_min_thread_budget]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:797-797`
- [[simulation.setup__rhs_flat_packet_min_items|_rhs_flat_packet_min_items]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:438-438`
- [[simulation.setup__rhs_flat_packet_overhead_min_samples|_rhs_flat_packet_overhead_min_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:454-454`
- [[simulation.setup__rhs_harmonics_batch_min_sats_per_worker|_rhs_harmonics_batch_min_sats_per_worker]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:819-819`
- [[simulation.setup__rhs_invsq_flat_min_sats|_rhs_invsq_flat_min_sats]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:963-963`
- [[simulation.setup__srp_ephemeris_cache_max_samples|_srp_ephemeris_cache_max_samples]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:201-201`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Clamping to 1 hides operator mistakes (a `-5` is accepted without warning). Values like `"1e3"` or `"4.0"` fail to parse as `Int` and throw. There is no upper bound, so absurdly large thresholds are accepted.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 23.
