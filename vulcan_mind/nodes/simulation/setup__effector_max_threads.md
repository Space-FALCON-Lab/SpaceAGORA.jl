---
id: simulation.setup__effector_max_threads
label: _effector_max_threads
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_max_threads
  lines:
  - 401
  - 401
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
  description: Return value of `_effector_max_threads`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _effector_max_threads

## Purpose
Caps how many threads the inner effector loop may use per RHS evaluation, limiting contention with the outer batch parallelism and with lock-guarded models such as GRAM.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_MAX_THREADS", 4)`, clamped to at least 1. The decision logic takes `min(available budget, this cap)` when computing the worker count for `_dynamic_effector_thread_decision`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_effector_max_threads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:857-857`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:402-402`
<!-- vulcan:connections:end -->

## Limitations
The cap is not validated against `Threads.nthreads()`, so a value above the pool is harmless but misleading in telemetry. A cap of 1 is equivalent to `:off` for the inner loop yet leaves the mode reporting `:auto`. Reads raw `ENV` rather than the engine override layer.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 401.
