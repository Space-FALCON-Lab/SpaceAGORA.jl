---
id: simulation.setup__effector_thread_threshold
label: _effector_thread_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_thread_threshold
  lines:
  - 397
  - 397
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
  description: Return value of `_effector_thread_threshold`.
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

# _effector_thread_threshold

## Purpose
Minimum number of effector work items (effectors × satellites in the flat view) required before the inner effector loop is considered for threading.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_EFFECTOR_THREAD_THRESHOLD", 2)`, clamped to at least 1. Feeds `_dynamic_effector_thread_decision`, where item counts below the threshold force a serial decision even in `:on` mode. The default of 2 is deliberately permissive so that cost-based checks, not the count, dominate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_effector_thread_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:856-856`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:398-398`
<!-- vulcan:connections:end -->

## Limitations
A threshold of 2 means the cost model (`_effector_work_ns_per_worker_threshold`) does almost all the work; if cost samples are unavailable early in a run the decision may flip-flop. No upper bound.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 397.
