---
id: parallel.context__spin_barrier_dispatch_bang
label: _spin_barrier_dispatch!
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _spin_barrier_dispatch!
  lines:
  - 208
  - 208
inputs:
- id: pool
  type: _SpinBarrierPool
  units: n/a
  required: true
  description: Positional argument `pool`.
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
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
  description: 'Return value of `_spin_barrier_dispatch!`; mutates `pool` in place.
    Returns `nothing`. Type parameters: `{F <: Function}`.'
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

# _spin_barrier_dispatch!

## Purpose
Executes one batch of `num_items` over a spin-barrier pool with `workers` logical workers, publishing the request, waking `workers - 1` spinning tasks, running the last worker slot on the coordinating thread itself, and then spin-waiting for completion and propagating errors.

## Design & Implementation
Computes `pool_workers = workers - 1`, reads `inner_scheduler_mode()` and `inner_dynamic_chunk_size()`, and allocates a fresh `Threads.Atomic{Int}(1)` cursor. Under `pool.run_lock` it stores a `NamedTuple` request in `pool.request[]`, then `atomic_add!`s `pool.worker_gen[w]` for `w in 1:pool_workers` to release them. The coordinator executes slot `workers` itself (dynamic chunks via the shared cursor or the stride `workers:workers:num_items`), catching any error into `coordinator_error`. It then spins with `GC.safepoint()` until `pool.done_count[] >= pool_workers`, subtracts `pool_workers` from the counter, rethrows the first non-`nothing` `pool.errors[w].ex` by worker index, and finally rethrows `coordinator_error`. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pool` | _SpinBarrierPool | n/a | yes | Positional argument `pool`. |
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `f` | F | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_spin_barrier_dispatch!`; mutates `pool` in place. Returns `nothing`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.thread_execution_threaded_foreach_worker_spin|threaded_foreach_worker_spin]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:173-173`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:243-243`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:243-243`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/context.jl:243-243`
- `callees` → [[parallel.env_config_inner_dynamic_chunk_size|inner_dynamic_chunk_size]] · `callers` · call · `src/parallel/policy/context.jl:219-219`
- `callees` → [[parallel.env_config_inner_scheduler_mode|inner_scheduler_mode]] · `callers` · call · `src/parallel/policy/context.jl:218-218`
<!-- vulcan:connections:end -->

## Limitations
Requires `workers - 1 <= pool.workers`; an oversized request indexes `worker_gen` out of bounds under `@inbounds`. The coordinator spins rather than yielding, so a long-running pool worker keeps the calling thread at 100 percent. Only the first pool-worker error is thrown; the backtrace is discarded by rethrowing `.ex`. The request `Ref` is shared, so a worker that lags a full generation could in principle read a newer request; the `run_lock` and barrier prevent this only if every worker completes each round.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 208.
