---
id: parallel.thread_execution_thread_worker_count
label: thread_worker_count
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: thread_worker_count
  lines:
  - 195
  - 195
inputs:
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: allotment
  type: Int
  units: n/a
  required: true
  description: Positional argument `allotment`.
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
  description: Return value of `thread_worker_count`.
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

# thread_worker_count

## Purpose
Public, exported wrapper around `_thread_worker_count` that lets modules outside `ParallelPolicy` (for example the GRAM isolated-pool density evaluator) compute the same worker count the foreach primitives will use, so they can size per-worker resources to match.

## Design & Implementation
An `@inline` function with signature `(num_items::Int, allotment::Int)::Int` that returns `_thread_worker_count(num_items, allotment)` unchanged. The private implementation applies `min(num_items, max(1, allotment), effective_inner_thread_budget())` and collapses to 1 when `Threads.nthreads() <= 1`. Keeping the public name separate from the private one allows the internal rule to change without breaking the exported API surface.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `thread_worker_count`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:775-775`
- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:927-927`
- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:995-995`
- [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:915-915`
- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:197-197`

**Downstream**

- `callees` → [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callers` · call · `src/parallel/policy/thread_execution.jl:196-196`
<!-- vulcan:connections:end -->

## Limitations
Because the answer depends on `effective_inner_thread_budget()` and the environment at call time, a caller that computes the count and then invokes a foreach primitive later may observe a mismatch if the budget changed in between. The function does not reserve threads; two callers each asking for the full count will both be told the same number. It cannot express the `callback_persistent_workers_enabled()` switch, which the persistent variants apply separately.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 195.
