---
id: core.runtime_types_rhsplanenvconfig
label: RhsPlanEnvConfig
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: RhsPlanEnvConfig
  lines:
  - 647
  - 647
inputs:
- id: execution_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `execution_mode`.
- id: profile_forces_serial
  type: Bool
  units: n/a
  required: true
  description: Field `profile_forces_serial`.
- id: batch_parallel_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `batch_parallel_mode`.
- id: batch_thread_threshold
  type: Int
  units: n/a
  required: true
  description: Field `batch_thread_threshold`.
- id: effector_parallel_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `effector_parallel_mode`.
- id: effector_thread_threshold
  type: Int
  units: n/a
  required: true
  description: Field `effector_thread_threshold`.
- id: effector_max_threads
  type: Int
  units: n/a
  required: true
  description: Field `effector_max_threads`.
- id: effector_allow_with_outer
  type: Bool
  units: n/a
  required: true
  description: Field `effector_allow_with_outer`.
- id: effector_heavy_only
  type: Bool
  units: n/a
  required: true
  description: Field `effector_heavy_only`.
- id: effector_cost_ns_per_item_default
  type: Float64
  units: n/a
  required: true
  description: Field `effector_cost_ns_per_item_default`.
- id: effector_cost_min_samples
  type: Int
  units: n/a
  required: true
  description: Field `effector_cost_min_samples`.
- id: effector_cost_ema_alpha
  type: Float64
  units: n/a
  required: true
  description: Field `effector_cost_ema_alpha`.
- id: effector_work_ns_per_worker_threshold
  type: Float64
  units: n/a
  required: true
  description: Field `effector_work_ns_per_worker_threshold`.
- id: effector_outer_work_scale
  type: Float64
  units: n/a
  required: true
  description: Field `effector_outer_work_scale`.
- id: flat_min_sats
  type: Int
  units: n/a
  required: true
  description: Field `flat_min_sats`.
- id: flat_min_effectors
  type: Int
  units: n/a
  required: true
  description: Field `flat_min_effectors`.
- id: flat_work_ns_threshold
  type: Float64
  units: n/a
  required: true
  description: Field `flat_work_ns_threshold`.
- id: flat_work_per_worker_ns_threshold
  type: Float64
  units: n/a
  required: true
  description: Field `flat_work_per_worker_ns_threshold`.
- id: flat_cost_heterogeneity_threshold
  type: Float64
  units: n/a
  required: true
  description: Field `flat_cost_heterogeneity_threshold`.
- id: flat_min_thread_budget
  type: Int
  units: n/a
  required: true
  description: Field `flat_min_thread_budget`.
- id: harmonics_batch_enabled
  type: Bool
  units: n/a
  required: true
  description: Field `harmonics_batch_enabled`.
- id: harmonics_batch_min_sats_per_worker
  type: Int
  units: n/a
  required: true
  description: Field `harmonics_batch_min_sats_per_worker`.
- id: harmonics_batch_spin_barrier
  type: Bool
  units: n/a
  required: true
  description: Field `harmonics_batch_spin_barrier`.
- id: harmonics_batch_allow_with_outer
  type: Bool
  units: n/a
  required: true
  description: Field `harmonics_batch_allow_with_outer`.
- id: rhs_effector_cost_min_samples
  type: Int
  units: n/a
  required: true
  description: Field `rhs_effector_cost_min_samples`.
- id: flat_packet_target_min_ns
  type: Float64
  units: n/a
  required: true
  description: Field `flat_packet_target_min_ns`.
- id: flat_packet_scheduler_mode
  type: Symbol
  units: n/a
  required: true
  description: Field `flat_packet_scheduler_mode`.
- id: flat_packet_min_items
  type: Int
  units: n/a
  required: true
  description: Field `flat_packet_min_items`.
- id: flat_packet_work_ns_threshold
  type: Float64
  units: n/a
  required: true
  description: Field `flat_packet_work_ns_threshold`.
- id: flat_packet_heterogeneity_threshold
  type: Float64
  units: n/a
  required: true
  description: Field `flat_packet_heterogeneity_threshold`.
- id: flat_packet_overhead_disable_ratio
  type: Float64
  units: n/a
  required: true
  description: Field `flat_packet_overhead_disable_ratio`.
- id: flat_packet_overhead_min_samples
  type: Int
  units: n/a
  required: true
  description: Field `flat_packet_overhead_min_samples`.
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
  type: RhsPlanEnvConfig
  units: n/a
  description: Constructed `RhsPlanEnvConfig`.
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

# RhsPlanEnvConfig

## Purpose
The run-scoped snapshot of the thirty-odd environment knobs that steer the per-RHS-call execution-plan routing in the simulation engine.

## Design & Implementation
An immutable struct grouping the execution mode and serial-profiling flag, batch and effector threading modes with thresholds and cost-model parameters, the flat-mode admission thresholds, harmonics batch settings, and the flat-packet scheduler parameters including its overhead-based self-disable ratio. Built once at setup and read as plain fields by `_rhs_execution_plan`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `execution_mode` | Symbol | n/a | yes | Field `execution_mode`. |
| in | `profile_forces_serial` | Bool | n/a | yes | Field `profile_forces_serial`. |
| in | `batch_parallel_mode` | Symbol | n/a | yes | Field `batch_parallel_mode`. |
| in | `batch_thread_threshold` | Int | n/a | yes | Field `batch_thread_threshold`. |
| in | `effector_parallel_mode` | Symbol | n/a | yes | Field `effector_parallel_mode`. |
| in | `effector_thread_threshold` | Int | n/a | yes | Field `effector_thread_threshold`. |
| in | `effector_max_threads` | Int | n/a | yes | Field `effector_max_threads`. |
| in | `effector_allow_with_outer` | Bool | n/a | yes | Field `effector_allow_with_outer`. |
| in | `effector_heavy_only` | Bool | n/a | yes | Field `effector_heavy_only`. |
| in | `effector_cost_ns_per_item_default` | Float64 | n/a | yes | Field `effector_cost_ns_per_item_default`. |
| in | `effector_cost_min_samples` | Int | n/a | yes | Field `effector_cost_min_samples`. |
| in | `effector_cost_ema_alpha` | Float64 | n/a | yes | Field `effector_cost_ema_alpha`. |
| in | `effector_work_ns_per_worker_threshold` | Float64 | n/a | yes | Field `effector_work_ns_per_worker_threshold`. |
| in | `effector_outer_work_scale` | Float64 | n/a | yes | Field `effector_outer_work_scale`. |
| in | `flat_min_sats` | Int | n/a | yes | Field `flat_min_sats`. |
| in | `flat_min_effectors` | Int | n/a | yes | Field `flat_min_effectors`. |
| in | `flat_work_ns_threshold` | Float64 | n/a | yes | Field `flat_work_ns_threshold`. |
| in | `flat_work_per_worker_ns_threshold` | Float64 | n/a | yes | Field `flat_work_per_worker_ns_threshold`. |
| in | `flat_cost_heterogeneity_threshold` | Float64 | n/a | yes | Field `flat_cost_heterogeneity_threshold`. |
| in | `flat_min_thread_budget` | Int | n/a | yes | Field `flat_min_thread_budget`. |
| in | `harmonics_batch_enabled` | Bool | n/a | yes | Field `harmonics_batch_enabled`. |
| in | `harmonics_batch_min_sats_per_worker` | Int | n/a | yes | Field `harmonics_batch_min_sats_per_worker`. |
| in | `harmonics_batch_spin_barrier` | Bool | n/a | yes | Field `harmonics_batch_spin_barrier`. |
| in | `harmonics_batch_allow_with_outer` | Bool | n/a | yes | Field `harmonics_batch_allow_with_outer`. |
| in | `rhs_effector_cost_min_samples` | Int | n/a | yes | Field `rhs_effector_cost_min_samples`. |
| in | `flat_packet_target_min_ns` | Float64 | n/a | yes | Field `flat_packet_target_min_ns`. |
| in | `flat_packet_scheduler_mode` | Symbol | n/a | yes | Field `flat_packet_scheduler_mode`. |
| in | `flat_packet_min_items` | Int | n/a | yes | Field `flat_packet_min_items`. |
| in | `flat_packet_work_ns_threshold` | Float64 | n/a | yes | Field `flat_packet_work_ns_threshold`. |
| in | `flat_packet_heterogeneity_threshold` | Float64 | n/a | yes | Field `flat_packet_heterogeneity_threshold`. |
| in | `flat_packet_overhead_disable_ratio` | Float64 | n/a | yes | Field `flat_packet_overhead_disable_ratio`. |
| in | `flat_packet_overhead_min_samples` | Int | n/a | yes | Field `flat_packet_overhead_min_samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RhsPlanEnvConfig | n/a | — | Constructed `RhsPlanEnvConfig`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:850-850`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Thirty-two fields with no grouping into sub-structs make the routing configuration hard to reason about as a whole; every knob is independent, so incompatible combinations are only discovered at plan time.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 647.
