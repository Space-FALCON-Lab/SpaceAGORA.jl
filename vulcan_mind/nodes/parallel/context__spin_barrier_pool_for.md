---
id: parallel.context__spin_barrier_pool_for
label: _spin_barrier_pool_for
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _spin_barrier_pool_for
  lines:
  - 283
  - 283
inputs:
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
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
  description: Return value of `_spin_barrier_pool_for`.
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

# _spin_barrier_pool_for

## Purpose
Returns the `_SpinBarrierPool` for a `source` symbol within the active policy scope, creating it on first use with `Threads.nthreads() - 1` workers so the coordinator always retains one thread.

## Design & Implementation
Builds the key with `_persistent_pool_key(source)` and, under `lock(_spin_barrier_lock)`, calls `get!(_spin_barrier_pools, key) do _create_spin_barrier_pool(max(1, Threads.nthreads() - 1)) end`. The closure runs only on a cache miss. The `return` inside the `do` block returns the pool from the function. Return type is annotated `_SpinBarrierPool`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _SpinBarrierPool | n/a | — | Return value of `_spin_barrier_pool_for`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.thread_execution_threaded_foreach_worker_spin|threaded_foreach_worker_spin]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:172-172`

**Downstream**

- `callees` → [[parallel.context__create_spin_barrier_pool|_create_spin_barrier_pool]] · `callers` · call · `src/parallel/policy/context.jl:289-289`
- `callees` → [[parallel.thread_execution__persistent_pool_key|_persistent_pool_key]] · `callers` · call · `src/parallel/policy/context.jl:284-284`
<!-- vulcan:connections:end -->

## Limitations
The pool size ignores `effective_inner_thread_budget()`, so under an outer-parallel split a spin pool still grabs all but one thread and its spinning workers will fight with sibling outer workers for cores. Pools live until the scope is destroyed. Because the key includes only scope and source, two callers with different `allotment` share the same pool and its fixed worker count.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 283.
