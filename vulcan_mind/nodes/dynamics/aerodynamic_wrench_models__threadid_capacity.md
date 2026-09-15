---
id: dynamics.aerodynamic_wrench_models__threadid_capacity
label: _threadid_capacity
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _threadid_capacity
  lines:
  - 20
  - 20
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
  description: Return value of `_threadid_capacity`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _threadid_capacity

## Purpose
Compatibility shim returning a scratch-buffer capacity large enough for either thread-id or worker-id indexing after the migration to stable worker ids.

## Design & Implementation
Returns `max(Threads.maxthreadid(), _multibody_max_threads())`. The comment notes it exists for legacy tests and callers that still size per-thread scratch by thread id.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_threadid_capacity`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_max_threads|_multibody_max_threads]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:23-23`
<!-- vulcan:connections:end -->

## Limitations
Nothing in this file calls it any more; it exists solely for external legacy callers. `Threads.maxthreadid()` can exceed `nthreads()` when interactive threads are enabled, over-allocating scratch.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 20.
