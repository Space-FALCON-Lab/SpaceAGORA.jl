---
id: parallel.thread_execution_threaded_collect_bang
label: threaded_collect!
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_collect!
  lines:
  - 252
  - 252
inputs:
- id: results
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `results`.
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
- id: f
  type: F
  units: n/a
  required: true
  description: Positional argument `f`.
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
  description: 'Return value of `threaded_collect!`; mutates `results` in place. Returns
    `nothing` or `results`. Type parameters: `{F <: Function}`.'
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

# threaded_collect!

## Purpose
Evaluates `f(idx)` for every item in parallel and stores each value at `results[idx]`, leaving all accumulation to the caller so that sums feeding the integrator are computed in deterministic index order independent of worker count or scheduling.

## Design & Implementation
Implemented as a closure over `threaded_foreach_worker(num_items, allotment) do _, idx ... end` that writes `@inbounds results[idx] = f(idx)` and returns `nothing`; the `worker_id` argument is ignored. The `results` vector is returned for convenience. A second method with `f` first supports `do`-block syntax. Because each index is written exactly once by exactly one worker, no synchronisation is required beyond the `@sync` inside the underlying primitive.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `results` | AbstractVector | n/a | yes | Positional argument `results`. |
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `f` | F | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `threaded_collect!`; mutates `results` in place. Returns `nothing` or `results`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.thread_execution_threaded_collect_persistent_bang|threaded_collect_persistent!]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:269-269`
- [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callees` → `callers` · feedback · `src/parallel/policy/thread_execution.jl:241-241`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_bang|_accumulate_dynamic_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:81-81`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:142-142`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:254-254`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:254-254`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:254-254`
- `callees` → [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callers` · call · `src/parallel/policy/thread_execution.jl:253-253`
<!-- vulcan:connections:end -->

## Limitations
`results` must already have length at least `num_items`; the `@inbounds` write makes an undersized vector undefined behaviour rather than a `BoundsError`. Element type of `results` must accept the return type of `f`, otherwise a conversion error is raised on the worker and propagates. The non-persistent primitive spawns fresh tasks each call, so for frequent small batches `threaded_collect_persistent!` is preferable.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 252.
