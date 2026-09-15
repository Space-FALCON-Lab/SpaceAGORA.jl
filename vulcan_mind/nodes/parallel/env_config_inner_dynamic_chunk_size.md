---
id: parallel.env_config_inner_dynamic_chunk_size
label: inner_dynamic_chunk_size
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: inner_dynamic_chunk_size
  lines:
  - 63
  - 63
inputs:
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
  description: Return value of `inner_dynamic_chunk_size`.
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

# inner_dynamic_chunk_size

## Purpose
Returns the number of work items a thread claims at a time when the inner scheduler is in `:dynamic` mode, trading load balance against per-claim synchronisation overhead.

## Design & Implementation
A one-line wrapper: `parse_thread_threshold_env("SPACEAGORA_PARALLEL_POLICY_INNER_DYNAMIC_CHUNK", 1)`. The default of 1 gives finest-grained balancing; the underlying parser clamps to at least 1 and throws `ArgumentError` on non-integer text. The value is meaningful only when `inner_scheduler_mode()` returns `:dynamic`, but the function does not check that.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `inner_dynamic_chunk_size`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.context__spin_barrier_dispatch_bang|_spin_barrier_dispatch!]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:219-219`
- [[parallel.thread_execution__threaded_foreach_persistent_bang|_threaded_foreach_persistent!]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:64-64`
- [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:210-210`
- [[parallel.thread_execution_threaded_reduce|threaded_reduce]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:300-300`
- [[parcore.thread_execution_threaded_foreach|threaded_foreach]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:12-12`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/parallel/policy/env_config.jl:64-64`
<!-- vulcan:connections:end -->

## Limitations
A chunk size larger than the item count degenerates to a single thread doing everything with no warning. There is no upper bound and no relation enforced between chunk size and thread count. Each call re-parses `ENV`.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 63.
