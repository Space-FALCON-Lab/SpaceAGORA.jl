---
id: parallel.env_config_inner_scheduler_mode
label: inner_scheduler_mode
kind: function
source:
  file: src/parallel/policy/env_config.jl
  symbol: inner_scheduler_mode
  lines:
  - 53
  - 53
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
  type: Symbol
  units: n/a
  description: Return value of `inner_scheduler_mode`.
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

# inner_scheduler_mode

## Purpose
Selects how inner-loop work (per-satellite or per-effector items within one RHS evaluation) is distributed across threads: static striding or dynamic chunked scheduling.

## Design & Implementation
Reads `SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER` with default `"static"`, lowercases and strips it. Both `"static"` and `"strided"` map to `:static`; `"dynamic"` maps to `:dynamic`. Any other value throws `ArgumentError` naming the two accepted modes. Returns a `Symbol` that the scheduler dispatches on; the dynamic mode additionally consults `inner_dynamic_chunk_size` for its chunk length.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `inner_scheduler_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.context__spin_barrier_dispatch_bang|_spin_barrier_dispatch!]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:218-218`
- [[parallel.thread_execution__threaded_foreach_persistent_bang|_threaded_foreach_persistent!]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:63-63`
- [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:208-208`
- [[parallel.thread_execution_threaded_reduce|threaded_reduce]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:298-298`
- [[parcore.thread_execution_threaded_foreach|threaded_foreach]] · `callees` → `callers` · call · `src/parallel/policy/thread_execution.jl:10-10`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only two modes exist; there is no `:auto` that would choose based on measured heterogeneity. Empty string values throw. The function re-reads `ENV` every call and is not snapshotted into `PolicyDecisionEnvConfig`, so changing the variable mid-run takes effect on the next call.

## Provenance
Mapped from `src/parallel/policy/env_config.jl` line 53.
