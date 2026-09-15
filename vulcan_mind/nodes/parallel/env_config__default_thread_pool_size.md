---
id: parallel.env_config__default_thread_pool_size
label: _default_thread_pool_size
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: _default_thread_pool_size
  lines:
  - 114
  - 114
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
  description: Return value of `_default_thread_pool_size`.
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

# _default_thread_pool_size

## Purpose
Reports the number of Julia threads available for inner parallel work, used as the ceiling when computing the effective thread budget.

## Design & Implementation
Wraps `Threads.nthreads()` in a `try`/`catch` whose fallback is also `Threads.nthreads()`, so both branches return the same value; the try block is effectively a no-op left over from an earlier implementation that may have queried a different pool. Returns `Int`, always at least 1 for a running Julia process.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_default_thread_pool_size`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/env_config.jl`
- [[parallel.env_config_effective_inner_thread_budget|effective_inner_thread_budget]] · `callees` → `callers` · call · `src/parallel/policy/env_config.jl:129-129`
- [[parallel.thread_execution__persistent_pool_for|_persistent_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:52-52`
- [[parallel.thread_execution__persistent_worker_pool_for|_persistent_worker_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:122-122`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `catch` branch is unreachable in practice because `Threads.nthreads()` does not throw, so the guard adds nothing. It reports the default thread pool only and ignores interactive-pool threads. Per-call evaluation is cheap but the result is nonetheless cached by `snapshot_policy_decision_env` via `effective_inner_thread_budget`.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 114.
