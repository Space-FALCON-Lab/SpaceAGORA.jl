---
id: dynamics.aerodynamic_wrench_models__multibody_max_threads
label: _multibody_max_threads
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _multibody_max_threads
  lines:
  - 16
  - 16
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
  description: Return value of `_multibody_max_threads`.
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

# _multibody_max_threads

## Purpose
Returns the cap on threads the multibody aerodynamic loop may use, read from `SPACEAGORA_MULTIBODY_MAX_THREADS` with default 4.

## Design & Implementation
Delegates to `ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_MAX_THREADS", 4)`. `_multibody_thread_decision` takes `min(policy.allotment, _multibody_max_threads())` as the final allotment, and `_threadid_capacity` uses it as a lower bound for scratch sizing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_multibody_max_threads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:51-51`
- [[dynamics.aerodynamic_wrench_models__threadid_capacity|_threadid_capacity]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:23-23`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:17-17`
<!-- vulcan:connections:end -->

## Limitations
The cap is not reconciled with `Threads.nthreads()`, so a value above the Julia thread count is silently clamped only by the policy layer. The parser shared with the threshold treats this as a threshold-style integer, so zero is not obviously rejected here.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 16.
