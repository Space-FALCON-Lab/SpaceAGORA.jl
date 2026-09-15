---
id: parcore.context_with_policy_context
label: with_policy_context
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: with_policy_context
  lines:
  - 332
  - 342
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ParallelPolicy namespace supplying the policy context record, the pool
    registries and their locks.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: scope
  type: Any
  units: n/a
  description: Return value of the wrapped function, produced with a fresh policy
    context installed and all worker pools created inside the scope torn down.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# with_policy_context

## Purpose
`with_policy_context` establishes an isolated policy scope for a block of work. Inside the block, telemetry, adaptive controller state and every persistent worker pool belong to a fresh `PolicyContext` stored in task-local storage, and when the block exits, all pools created under that scope are shut down. It is what makes a nested campaign or a test case independent of whatever ran before it.

## Model & Assumptions
Scope identity is a `UInt` derived from the context object, and pool registries are keyed by the pair of scope id and source symbol. That keying is the isolation mechanism: two scopes asking for the same source get different pools. The context is installed with `Base.task_local_storage`, so a spawned task inherits the storage of the task that spawned it, and code running outside any scope falls back to `_global_policy_context`.

## Design & Implementation
The function creates a `PolicyContext`, computes its scope id, and calls `Base.task_local_storage` with the key `:spaceagora_parallel_policy_context`. The wrapped call is inside a `try`/`finally` whose cleanup calls `_destroy_persistent_foreach_scope!`, so pools are released even when the body throws. That destroy routine collects stale keys from `_persistent_foreach_pools` and `_persistent_foreach_worker_pools` under `_persistent_foreach_lock`, and from `_spin_barrier_pools` under `_spin_barrier_lock`, deletes them from the registries, and only then shuts each pool down outside the locks — the ordering avoids holding a registry lock while waiting on worker tasks. The same file defines the channel-based persistent foreach pools, the spin-barrier pool with its dispatch and shutdown paths, and the per-worker loops those pools run.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ParallelPolicy namespace supplying the policy context record, the pool registries and their locks. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `scope` | Any | n/a | — | Return value of the wrapped function, produced with a fresh policy context installed and all worker pools created inside the scope torn down. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:152-152`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:337-337`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/context.jl:337-337`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/context.jl:337-337`
- `callees` → [[parallel.context__destroy_persistent_foreach_scope_bang|_destroy_persistent_foreach_scope!]] · `callers` · call · `src/parallel/policy/context.jl:339-339`
- `callees` → [[parallel.context__policy_scope_id|_policy_scope_id]] · `callers` · call · `src/parallel/policy/context.jl:334-334`
- `callees` → [[parallel.types_policycontext|PolicyContext]] · `callers` · call · `src/parallel/policy/context.jl:333-333`
<!-- vulcan:connections:end -->

## Limitations
Isolation is per task tree, not per thread: work handed to a pool that predates the scope keeps the older context. Because cleanup happens only at scope exit, a long-lived scope that touches many sources accumulates pools and their worker tasks for its whole lifetime. Shutting down a spin-barrier pool requires its workers to observe the stop flag, so teardown of a busy pool is not instantaneous.

## Provenance
Mapped from `src/parallel/policy/context.jl:332-342`.
