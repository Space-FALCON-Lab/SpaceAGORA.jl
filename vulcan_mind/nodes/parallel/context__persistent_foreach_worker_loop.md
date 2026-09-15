---
id: parallel.context__persistent_foreach_worker_loop
label: _persistent_foreach_worker_loop
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _persistent_foreach_worker_loop
  lines:
  - 21
  - 21
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
  description: Return value of `_persistent_foreach_worker_loop`.
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

# _persistent_foreach_worker_loop

## Purpose
Body of each long-lived task in an index-only persistent pool: blocks on its request channel, executes its share of a dispatched batch by calling `f(idx)`, reports completion or a captured exception on the shared done channel, and exits on `:stop`.

## Design & Implementation
An infinite `while true` loop performing `request = take!(request_channel)`; the sentinel `:stop` returns `nothing`. Otherwise the request `NamedTuple` fields `num_items`, `active_workers`, `scheduler`, `chunk`, and `f` are unpacked inside a `try`. In `:dynamic` mode the shared `request.next_index` atomic hands out `chunk`-sized index ranges via `atomic_add!` until `start_idx > num_items`; in static mode the worker iterates `worker_id:active_workers:num_items`. Any exception is wrapped as `Base.CapturedException(err, catch_backtrace())`. The loop always `put!`s either `nothing` or the captured exception into `done_channel`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `worker_id` | Int | n/a | yes | Positional argument `worker_id`. |
| in | `request_channel` | Channel{Any} | n/a | yes | Positional argument `request_channel`. |
| in | `done_channel` | Channel{Any} | n/a | yes | Positional argument `done_channel`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_persistent_foreach_worker_loop`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/context.jl`
- [[parallel.context__create_persistent_foreach_pool|_create_persistent_foreach_pool]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:130-130`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:46-46`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:46-46`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/context.jl:46-46`
<!-- vulcan:connections:end -->

## Limitations
Because `request` is typed `Any`, field access and the call to `f` are dynamically dispatched, adding overhead per item that is negligible for heavy bodies but significant for very cheap ones. A worker with `worker_id > active_workers` in static mode still receives a request and runs an empty range; the dispatcher only sends to the first `workers` channels, so idle workers stay blocked. If `f` never returns the pool is permanently wedged.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 21.
