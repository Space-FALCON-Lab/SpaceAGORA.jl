---
id: parallel.context__shutdown_spin_barrier_pool_bang
label: _shutdown_spin_barrier_pool!
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _shutdown_spin_barrier_pool!
  lines:
  - 274
  - 274
inputs:
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
  description: Return value of `_shutdown_spin_barrier_pool!`; mutates `pool` in place.
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

# _shutdown_spin_barrier_pool!

## Purpose
Stops all busy-polling workers of a `_SpinBarrierPool` by raising the shared `stop` flag and bumping every worker's generation counter so each spinner wakes, observes the flag, and returns.

## Design & Implementation
Sets `pool.stop[] = true`, then loops over `pool.worker_gen` and applies `Threads.atomic_add!(gen, 1)` to each. In `_spin_barrier_worker_loop_w` the inner poll loop checks `pool.stop[]` on every iteration and again immediately after a generation change, so the bump guarantees prompt exit even for a worker that was just about to read a request. Returns `nothing`. No lock is taken.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pool` | _SpinBarrierPool | n/a | yes | Positional argument `pool`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_shutdown_spin_barrier_pool!`; mutates `pool` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.context__destroy_persistent_foreach_scope_bang|_destroy_persistent_foreach_scope!]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:327-327`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `run_lock` is not acquired, calling this during an active `_spin_barrier_dispatch!` lets workers exit mid-batch without incrementing `done_count`, so the coordinator's spin-wait `while pool.done_count[] < pool_workers` never terminates. The scope-destruction path only calls this after the scope's user code has finished, avoiding that scenario. The function does not join the worker tasks.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 274.
