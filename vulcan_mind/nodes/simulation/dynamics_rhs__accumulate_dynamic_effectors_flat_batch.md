---
id: simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang
label: _accumulate_dynamic_effectors_flat_batch!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_dynamic_effectors_flat_batch!
  lines:
  - 983
  - 983
inputs:
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: plan
  type: Any
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: partition
  type: Union{Nothing, Symbol}
  units: n/a
  required: false
  description: Keyword argument `partition` (default `nothing`).
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
  type: Nothing
  units: n/a
  description: Return value of `_accumulate_dynamic_effectors_flat_batch!`; mutates
    `sc_state` in place. Returns `_accumulate_harmonics_flat_batch!(sc_state, p, t,
    dynamic_effectors[1], plan)` or `nothing`.
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

# _accumulate_dynamic_effectors_flat_batch!

## Purpose
The flat constellation effector queue: evaluates every satellite-by-effector work item across a worker pool, with batchable effectors handled by vectorised kernels and the rest by a dynamic queue.

## Design & Implementation
Computes worker count from the plan, ensures scratch buffers, and short-circuits to the harmonics batch kernel for a lone harmonics model. It prefills flat state samples, runs batchable effectors through `_accumulate_batchable_effector_flat!`, builds the execution plan and work items, optionally packets them, and dispatches through the persistent worker pool with per-worker partial buffers that are reduced into `totals`. Packet timings feed the cost model.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `plan` | Any | n/a | yes | Positional argument `plan`. |
| in | `partition` | Union{Nothing, Symbol} | n/a | no | Keyword argument `partition` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_accumulate_dynamic_effectors_flat_batch!`; mutates `sc_state` in place. Returns `_accumulate_harmonics_flat_batch!(sc_state, p, t, dynamic_effectors[1], plan)` or `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1332-1332`

**Downstream**

- `callees` → [[parallel.thread_execution_thread_worker_count|thread_worker_count]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:995-995`
- `callees` → [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1100-1100`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1169-1169`
- `callees` → [[simulation.dynamics_rhs__accumulate_batchable_effector_flat_bang|_accumulate_batchable_effector_flat!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1045-1045`
- `callees` → [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1010-1010`
- `callees` → [[simulation.dynamics_rhs__build_constellation_execution_plan_bang|_build_constellation_execution_plan!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1067-1067`
- `callees` → [[simulation.dynamics_rhs__constellation_node_eff_idx|_constellation_node_eff_idx]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1105-1105`
- `callees` → [[simulation.dynamics_rhs__constellation_node_sat_idx|_constellation_node_sat_idx]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1104-1104`
- `callees` → [[simulation.dynamics_rhs__count_flat_queue_only_effectors|_count_flat_queue_only_effectors]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1002-1002`
- `callees` → [[simulation.dynamics_rhs__count_non_batchable_effectors|_count_non_batchable_effectors]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1047-1047`
- `callees` → [[simulation.dynamics_rhs__ensure_rhs_flat_effector_scratch_bang|_ensure_rhs_flat_effector_scratch!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1003-1003`
- `callees` → [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1112-1112`
- `callees` → [[simulation.dynamics_rhs__has_any_batchable_effector|_has_any_batchable_effector]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1024-1024`
- `callees` → [[simulation.dynamics_rhs__has_any_harmonics_effector|_has_any_harmonics_effector]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1056-1056`
- `callees` → [[simulation.dynamics_rhs__partition_needs_state_sample|_partition_needs_state_sample]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1026-1026`
- `callees` → [[simulation.dynamics_rhs__partition_selected_count|_partition_selected_count]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1004-1004`
- `callees` → [[simulation.dynamics_rhs__prefill_rhs_flat_state_samples_bang|_prefill_rhs_flat_state_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1029-1029`
- `callees` → [[simulation.dynamics_rhs__prepare_rhs_flat_work_packets_bang|_prepare_rhs_flat_work_packets!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1084-1084`
- `callees` → [[simulation.dynamics_rhs__rhs_flat_state_sample_from_buffers|_rhs_flat_state_sample_from_buffers]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1110-1110`
- `callees` → [[simulation.dynamics_rhs__rhs_flat_use_packet_scheduler|_rhs_flat_use_packet_scheduler]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1078-1078`
- `callees` → [[simulation.dynamics_rhs__update_rhs_flat_packet_cost_model_bang|_update_rhs_flat_packet_cost_model!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1178-1178`
- `callees` → [[simulation.dynamics_rhs__update_rhs_flat_packet_overhead_model_bang|_update_rhs_flat_packet_overhead_model!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1189-1189`
- `callees` → [[simulation.dynamics_rhs__with_packet_scheduler|_with_packet_scheduler]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1079-1079`
- `callees` → [[simulation.effector_sampling__wrench_method_available|_wrench_method_available]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1023-1023`
- `callees` → [[simulation.setup__ensure_rhs_effector_cost_model_bang|_ensure_rhs_effector_cost_model!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1020-1020`
- `callees` → [[simulation.setup__update_effector_cost_model_bang|_update_effector_cost_model!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1163-1163`
<!-- vulcan:connections:end -->

## Limitations
Around 340 lines coordinating four scheduling strategies; the per-worker partials use a fixed six-row layout that must match `_flat_totals_force_torque`.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 983.
