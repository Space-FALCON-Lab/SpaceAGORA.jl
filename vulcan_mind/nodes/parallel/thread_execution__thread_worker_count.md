---
id: parallel.thread_execution__thread_worker_count
label: _thread_worker_count
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: _thread_worker_count
  lines:
  - 185
  - 185
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
  description: Return value of `_thread_worker_count`.
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

# _thread_worker_count

## Purpose
Central rule for how many worker threads an inner parallel loop should use given the item count, the caller's allotment hint, and the process-wide inner thread budget, returning 1 whenever parallelism is pointless or disabled.

## Design & Implementation
Returns 1 for `num_items <= 0`. Otherwise `budget = effective_inner_thread_budget()` (which honours `SPACEAGORA_INNER_THREAD_BUDGET` and outer-parallel splits) and `workers = min(num_items, max(1, allotment), budget)`. If that is 1 or fewer, or `Threads.nthreads() <= 1`, it returns 1; otherwise `workers`. The public wrapper `thread_worker_count` simply forwards to this `@inline` function so that external modules such as `ModelSelection` can size their own pools consistently.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_thread_worker_count`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.thread_execution_thread_worker_count|thread_worker_count]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:196-196`
- [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:101-101`
- [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:201-201`
- [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:134-134`
- [[parallel.thread_execution_threaded_foreach_worker_spin|threaded_foreach_worker_spin]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:165-165`
- [[parallel.thread_execution_threaded_reduce|threaded_reduce]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:284-284`
- [[parcore.thread_execution_threaded_foreach|threaded_foreach]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:3-3`

**Downstream**

- `callees` → [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callers` · call · `src/parallel/policy/thread_execution.jl:187-187`
<!-- vulcan:connections:end -->

## Limitations
`effective_inner_thread_budget()` may parse environment variables on each call unless cached, so this is not free in tight loops. A negative or zero `allotment` is clamped to 1, silently forcing serial execution rather than signalling misuse. The function does not consider current load or whether other pools are already busy, so nested callers can still oversubscribe if they each independently read the full budget.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 185.
