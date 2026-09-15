---
id: simulation.setup__snapshot_rhs_plan_env_config
label: _snapshot_rhs_plan_env_config
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _snapshot_rhs_plan_env_config
  lines:
  - 849
  - 849
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
  type: SimulationModel.RhsPlanEnvConfig
  units: n/a
  description: Return value of `_snapshot_rhs_plan_env_config`.
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

# _snapshot_rhs_plan_env_config

## Purpose
Resolves every environment knob the RHS execution planner consults into one typed `RhsPlanEnvConfig`, built once per run.

## Design & Implementation
Calls the thirty-two individual `_rhs_*`, `_effector_*` and `_profile_*` environment readers in the struct's field order and constructs the snapshot. Built inside any active `SimulationEngineConfig` override scope so overrides are captured.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.RhsPlanEnvConfig | n/a | — | Return value of `_snapshot_rhs_plan_env_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__harmonics_batch_allow_with_outer|_harmonics_batch_allow_with_outer]] · `callees` → `callers` · feedback · `src/simulation/engine/setup.jl:841-841`
- [[simulation.setup__initialize_runtime_env_config_bang|_initialize_runtime_env_config!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:914-914`
- [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:894-894`

**Downstream**

- `callees` → [[core.runtime_types_rhsplanenvconfig|RhsPlanEnvConfig]] · `callers` · call · `src/simulation/engine/setup.jl:850-850`
- `callees` → [[simulation.setup__effector_allow_with_outer|_effector_allow_with_outer]] · `callers` · call · `src/simulation/engine/setup.jl:858-858`
- `callees` → [[simulation.setup__effector_cost_ema_alpha|_effector_cost_ema_alpha]] · `callers` · call · `src/simulation/engine/setup.jl:862-862`
- `callees` → [[simulation.setup__effector_cost_min_samples|_effector_cost_min_samples]] · `callers` · call · `src/simulation/engine/setup.jl:861-861`
- `callees` → [[simulation.setup__effector_cost_ns_per_item_default|_effector_cost_ns_per_item_default]] · `callers` · call · `src/simulation/engine/setup.jl:860-860`
- `callees` → [[simulation.setup__effector_heavy_only|_effector_heavy_only]] · `callers` · call · `src/simulation/engine/setup.jl:859-859`
- `callees` → [[simulation.setup__effector_max_threads|_effector_max_threads]] · `callers` · call · `src/simulation/engine/setup.jl:857-857`
- `callees` → [[simulation.setup__effector_outer_work_scale|_effector_outer_work_scale]] · `callers` · call · `src/simulation/engine/setup.jl:864-864`
- `callees` → [[simulation.setup__effector_parallel_mode|_effector_parallel_mode]] · `callers` · call · `src/simulation/engine/setup.jl:855-855`
- `callees` → [[simulation.setup__effector_thread_threshold|_effector_thread_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:856-856`
- `callees` → [[simulation.setup__effector_work_ns_per_worker_threshold|_effector_work_ns_per_worker_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:863-863`
- `callees` → [[simulation.setup__harmonics_batch_allow_with_outer|_harmonics_batch_allow_with_outer]] · `callers` · call · `src/simulation/engine/setup.jl:874-874`
- `callees` → [[simulation.setup__profile_forces_serial_rhs|_profile_forces_serial_rhs]] · `callers` · call · `src/simulation/engine/setup.jl:852-852`
- `callees` → [[simulation.setup__rhs_batch_parallel_mode|_rhs_batch_parallel_mode]] · `callers` · call · `src/simulation/engine/setup.jl:853-853`
- `callees` → [[simulation.setup__rhs_batch_thread_threshold|_rhs_batch_thread_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:854-854`
- `callees` → [[simulation.setup__rhs_effector_cost_min_samples|_rhs_effector_cost_min_samples]] · `callers` · call · `src/simulation/engine/setup.jl:875-875`
- `callees` → [[simulation.setup__rhs_execution_mode_env|_rhs_execution_mode_env]] · `callers` · call · `src/simulation/engine/setup.jl:851-851`
- `callees` → [[simulation.setup__rhs_flat_cost_heterogeneity_threshold|_rhs_flat_cost_heterogeneity_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:869-869`
- `callees` → [[simulation.setup__rhs_flat_min_effectors|_rhs_flat_min_effectors]] · `callers` · call · `src/simulation/engine/setup.jl:866-866`
- `callees` → [[simulation.setup__rhs_flat_min_sats|_rhs_flat_min_sats]] · `callers` · call · `src/simulation/engine/setup.jl:865-865`
- `callees` → [[simulation.setup__rhs_flat_min_thread_budget|_rhs_flat_min_thread_budget]] · `callers` · call · `src/simulation/engine/setup.jl:870-870`
- `callees` → [[simulation.setup__rhs_flat_packet_heterogeneity_threshold|_rhs_flat_packet_heterogeneity_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:880-880`
- `callees` → [[simulation.setup__rhs_flat_packet_min_items|_rhs_flat_packet_min_items]] · `callers` · call · `src/simulation/engine/setup.jl:878-878`
- `callees` → [[simulation.setup__rhs_flat_packet_overhead_disable_ratio|_rhs_flat_packet_overhead_disable_ratio]] · `callers` · call · `src/simulation/engine/setup.jl:881-881`
- `callees` → [[simulation.setup__rhs_flat_packet_overhead_min_samples|_rhs_flat_packet_overhead_min_samples]] · `callers` · call · `src/simulation/engine/setup.jl:882-882`
- `callees` → [[simulation.setup__rhs_flat_packet_scheduler_mode|_rhs_flat_packet_scheduler_mode]] · `callers` · call · `src/simulation/engine/setup.jl:877-877`
- `callees` → [[simulation.setup__rhs_flat_packet_target_min_ns|_rhs_flat_packet_target_min_ns]] · `callers` · call · `src/simulation/engine/setup.jl:876-876`
- `callees` → [[simulation.setup__rhs_flat_packet_work_ns_threshold|_rhs_flat_packet_work_ns_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:879-879`
- `callees` → [[simulation.setup__rhs_flat_work_ns_threshold|_rhs_flat_work_ns_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:867-867`
- `callees` → [[simulation.setup__rhs_flat_work_per_worker_ns_threshold|_rhs_flat_work_per_worker_ns_threshold]] · `callers` · call · `src/simulation/engine/setup.jl:868-868`
- `callees` → [[simulation.setup__rhs_harmonics_batch_enabled|_rhs_harmonics_batch_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:871-871`
- `callees` → [[simulation.setup__rhs_harmonics_batch_min_sats_per_worker|_rhs_harmonics_batch_min_sats_per_worker]] · `callers` · call · `src/simulation/engine/setup.jl:872-872`
<!-- vulcan:connections:end -->

## Limitations
Positional construction means adding a field to the struct requires editing this call in the matching position, with no keyword safety.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 849.
