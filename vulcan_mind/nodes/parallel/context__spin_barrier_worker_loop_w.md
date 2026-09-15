---
id: parallel.context__spin_barrier_worker_loop_w
label: _spin_barrier_worker_loop_w
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _spin_barrier_worker_loop_w
  lines:
  - 152
  - 152
inputs:
- id: worker_id
  type: Int
  units: n/a
  required: true
  description: Positional argument `worker_id`.
- id: pool
  type: _SpinBarrierPool
  units: n/a
  required: true
  description: Positional argument `pool`.
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
  description: Return value of `_spin_barrier_worker_loop_w`.
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

# _spin_barrier_worker_loop_w

## Purpose
Body of each spin-barrier pool task: busy-polls its private generation counter until the coordinator bumps it, executes its share of the published request with `f(worker_id, idx)`, records success or a captured exception in its error slot, and increments the shared completion counter.

## Design & Implementation
Keeps a local `my_gen` starting at 0. The inner `while pool.worker_gen[worker_id][] == my_gen` loop checks `pool.stop[]` and calls `GC.safepoint()` on each iteration so garbage collection can proceed while spinning. On wake it increments `my_gen`, re-checks `stop`, reads `pool.request[]`, and inside a `try` processes either dynamic chunks from `request.next_index` or the static stride `worker_id:active_workers:num_items`, calling `f(worker_id, idx)`. Errors are wrapped in `Base.CapturedException`. It then writes `pool.errors[worker_id] = captured` and performs `Threads.atomic_add!(pool.done_count, 1)`, whose release semantics publish the slot write to the coordinator.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `worker_id` | Int | n/a | yes | Positional argument `worker_id`. |
| in | `pool` | _SpinBarrierPool | n/a | yes | Positional argument `pool`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_spin_barrier_worker_loop_w`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.context__create_spin_barrier_pool|_create_spin_barrier_pool]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:203-203`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:178-178`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:178-178`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/context.jl:178-178`
<!-- vulcan:connections:end -->

## Limitations
The task never yields, so it monopolises its Julia thread for the pool's lifetime; on a machine with fewer physical cores than `nthreads()` this causes heavy contention. A worker that misses a generation increment (two bumps before it wakes) would skip one batch and desynchronise `done_count`, hanging the coordinator; the `run_lock` prevents concurrent dispatches but not this theoretical lag under extreme scheduling delay. The `request` read is a plain `Ref` load ordered only by the atomic generation bump.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 152.
