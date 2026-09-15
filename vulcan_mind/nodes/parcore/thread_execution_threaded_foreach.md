---
id: parcore.thread_execution_threaded_foreach
label: threaded_foreach
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_foreach
  lines:
  - 1
  - 35
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelPolicy namespace supplying the scheduler mode, chunk size and
    worker-count helpers.
- id: allotment
  type: Int
  units: workers
  required: true
  description: Worker count decided by thread_policy_decision; a value of one or a
    single-threaded process forces the serial path.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: completion
  type: Nothing
  units: n/a
  description: 'Side effect only: the supplied closure has been applied exactly once
    to every index in 1:num_items before the call returns.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# threaded_foreach

## Purpose
`threaded_foreach` is the basic parallel loop primitive. It applies a closure to each index of a range using the worker count the policy allotted, and it is the execution half of the decide-execute-observe split that the policy module implements.

## Model & Assumptions
The contract is that the closure is called exactly once per index and that the call returns only after every index has been processed, which the `Threads.@sync` block enforces. Iterations must be independent: the primitive provides no ordering guarantee and no reduction, and any shared accumulation is the caller's responsibility. Non-positive item counts return immediately.

## Design & Implementation
The function first computes the effective worker count with `_thread_worker_count`, and falls back to a plain `@inbounds` serial loop when that is one or fewer or when the process has a single thread — this avoids paying task-spawn cost for work that cannot benefit. With more workers it selects between two schedulers read from the environment. The static scheduler spawns one task per worker and gives worker `w` the strided index set `w:workers:num_items`, which balances well when per-item cost is uniform and needs no coordination at all. The dynamic scheduler instead shares a `Threads.Atomic{Int}` cursor; each worker repeatedly claims a chunk with `atomic_add!`, stops when the claimed start index passes the item count, and otherwise processes the chunk, which handles ragged per-item cost at the price of one atomic operation per chunk. A second method takes the closure first so callers can use do-block syntax. The rest of the file layers the persistent, per-worker persistent, spin-barrier, collect and reduce variants on the same pattern.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelPolicy namespace supplying the scheduler mode, chunk size and worker-count helpers. |
| in | `allotment` | Int | workers | yes | Worker count decided by thread_policy_decision; a value of one or a single-threaded process forces the serial path. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `completion` | Nothing | n/a | — | Side effect only: the supplied closure has been applied exactly once to every index in 1:num_items before the call returns. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:103-103`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1341-1341`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:6-6`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:6-6`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:6-6`
- `callees` → [[parallel.env_config_inner_dynamic_chunk_size|inner_dynamic_chunk_size]] · `callers` · call · `src/parallel/policy/thread_execution.jl:12-12`
- `callees` → [[parallel.env_config_inner_scheduler_mode|inner_scheduler_mode]] · `callers` · call · `src/parallel/policy/thread_execution.jl:10-10`
- `callees` → [[parallel.thread_execution__thread_worker_count|_thread_worker_count]] · `callers` · call · `src/parallel/policy/thread_execution.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
Strided static scheduling is the default and is a poor fit when item cost correlates with index, since one worker inherits all the expensive entries. The dynamic path's chunk size is a single environment-wide value, so it cannot adapt per call site. Exceptions thrown inside a spawned task surface as a composite exception from the sync block, which obscures which index failed; the primitive does not attempt to report the failing index.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl:1-35`.
