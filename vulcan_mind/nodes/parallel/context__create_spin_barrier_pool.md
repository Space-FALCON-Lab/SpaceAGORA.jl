---
id: parallel.context__create_spin_barrier_pool
label: _create_spin_barrier_pool
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _create_spin_barrier_pool
  lines:
  - 197
  - 197
inputs:
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
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
  type: _SpinBarrierPool
  units: n/a
  description: Return value of `_create_spin_barrier_pool`.
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

# _create_spin_barrier_pool

## Purpose
Creates a `_SpinBarrierPool` and launches its busy-polling worker tasks, reserving one Julia thread for the coordinating task so that spinning workers can never starve the dispatcher.

## Design & Implementation
Clamps `workers = max(1, min(workers, Threads.nthreads() - 1))`, constructs `_SpinBarrierPool(workers)` (which allocates per-worker generation counters, error slots, a shared request `Ref`, a `done_count` atomic, and a `stop` flag), then `Threads.@spawn`s `_spin_barrier_worker_loop_w(worker_id, pool)` for each worker. The spawned tasks immediately enter their spin loop waiting for a generation bump. Returns the pool.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _SpinBarrierPool | n/a | — | Return value of `_create_spin_barrier_pool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.context__spin_barrier_pool_for|_spin_barrier_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:289-289`

**Downstream**

- `callees` → [[parallel.context__spin_barrier_worker_loop_w|_spin_barrier_worker_loop_w]] · `callers` · call · `src/parallel/policy/context.jl:203-203`
- `callees` → [[parallel.types__spinbarrierpool|_SpinBarrierPool]] · `callers` · call · `src/parallel/policy/context.jl:201-201`
<!-- vulcan:connections:end -->

## Limitations
Because workers spin without yielding, every pool worker permanently occupies a Julia thread; creating this pool on a machine where other threaded work must run degrades that work severely. With `Threads.nthreads() == 1` the clamp still yields 1 worker, which would deadlock the coordinator; `_spin_barrier_dispatch!` avoids this only because `workers - 1 == 0` pool tasks are signalled. Tasks are only released by `_shutdown_spin_barrier_pool!`.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 197.
