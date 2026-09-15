---
id: parallel.thread_execution__persistent_pool_key
label: _persistent_pool_key
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: _persistent_pool_key
  lines:
  - 44
  - 44
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
  type: Tuple{UInt,
  units: n/a
  description: Return value of `_persistent_pool_key`.
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

# _persistent_pool_key

## Purpose
Constructs the dictionary key under which persistent worker pools are stored, combining the current policy scope identifier with the caller-supplied `source` symbol so that pools are isolated per scope and per call site.

## Design & Implementation
An `@inline` function returning `(_active_policy_scope_id(), source)::Tuple{UInt, Symbol}`. `_active_policy_scope_id()` reads the task-local or global policy context so that nested or concurrent policy scopes (for example one per outer Monte Carlo worker) get independent pools even when they use the same `source` name. Both tuple elements are `isbits`, giving a cheap, allocation-free hash key.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{UInt, | n/a | — | Return value of `_persistent_pool_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.context__spin_barrier_pool_for|_spin_barrier_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:284-284`
- [[parallel.thread_execution__persistent_pool_for|_persistent_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:49-49`
- [[parallel.thread_execution__persistent_worker_pool_for|_persistent_worker_pool_for]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:119-119`

**Downstream**

- `callees` → [[parallel.context__active_policy_scope_id|_active_policy_scope_id]] · `callers` · call · `src/parallel/policy/thread_execution.jl:45-45`
<!-- vulcan:connections:end -->

## Limitations
Because the scope id is looked up at call time, calling from a task that has inherited a different policy context than expected silently selects a different pool. `Symbol` sources generated dynamically (via `Symbol(string)`) are interned forever, so key cardinality is bounded only by caller discipline. The key does not include the worker count, so a pool created under one thread budget is reused under another.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 44.
