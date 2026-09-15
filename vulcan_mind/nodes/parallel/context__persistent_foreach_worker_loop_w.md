---
id: parallel.context__persistent_foreach_worker_loop_w
label: _persistent_foreach_worker_loop_w
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _persistent_foreach_worker_loop_w
  lines:
  - 61
  - 61
inputs:
- id: worker_id
  type: Int
  units: n/a
  required: true
  description: Positional argument `worker_id`.
- id: request_channel
  type: Channel{Any}
  units: n/a
  required: true
  description: Positional argument `request_channel`.
- id: done_channel
  type: Channel{Any}
  units: n/a
  required: true
  description: Positional argument `done_channel`.
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
  description: Return value of `_persistent_foreach_worker_loop_w`.
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

# _persistent_foreach_worker_loop_w

## Purpose
Worker-aware counterpart of `_persistent_foreach_worker_loop`: identical dispatch, scheduling, and error-capture behaviour but invokes the user closure as `f(worker_id, idx)` so per-worker scratch state can be selected without locking.

## Design & Implementation
Loops on `take!(request_channel)`, exiting on the `:stop` sentinel. Inside a `try`, it reads `num_items`, `active_workers`, `scheduler`, `chunk`, and `f` from the request. In `:dynamic` mode it claims index ranges of length `chunk` from `request.next_index` with `Threads.atomic_add!` until exhaustion; otherwise it strides `worker_id:active_workers:num_items`. Each item is processed with `f(worker_id, idx)`. Exceptions become `Base.CapturedException` values. After every batch, `put!(done_channel, captured)` signals the coordinator, where `captured` is `nothing` on success.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `worker_id` | Int | n/a | yes | Positional argument `worker_id`. |
| in | `request_channel` | Channel{Any} | n/a | yes | Positional argument `request_channel`. |
| in | `done_channel` | Channel{Any} | n/a | yes | Positional argument `done_channel`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_persistent_foreach_worker_loop_w`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.context__create_persistent_foreach_worker_pool|_create_persistent_foreach_worker_pool]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:111-111`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:86-86`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:86-86`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/context.jl:86-86`
<!-- vulcan:connections:end -->

## Limitations
Shares all limitations of the index-only loop: `Any`-typed requests, dynamic dispatch of `f`, permanent blocking if `f` hangs, and only one captured exception per worker per batch. The `worker_id` is the pool slot, fixed at creation, so it is stable across batches but bears no relation to `Threads.threadid()`; callers must not use it to index thread-local arrays sized by `nthreads()`.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 61.
