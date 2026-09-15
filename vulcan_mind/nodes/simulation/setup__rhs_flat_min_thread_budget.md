---
id: simulation.setup__rhs_flat_min_thread_budget
label: _rhs_flat_min_thread_budget
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_min_thread_budget
  lines:
  - 796
  - 796
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
  description: Return value of `_rhs_flat_min_thread_budget`.
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

# _rhs_flat_min_thread_budget

## Purpose
Reads the smallest inner thread budget under which the flat effector queue is never used, because with too few workers its dynamic scheduling overhead outweighs the balance it buys.

## Design & Implementation
Parses `SPACEAGORA_EFFECTOR_FLAT_MIN_THREAD_BUDGET` through `parse_thread_threshold_env`, defaulting to the policy's `auto_thread_min_budget()` so the two thresholds stay aligned. Declared `@inline`. Captured into the RHS plan snapshot.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_flat_min_thread_budget`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:870-870`

**Downstream**

- `callees` → [[parallel.env_config_auto_thread_min_budget|auto_thread_min_budget]] · `callers` · call · `src/simulation/engine/setup.jl:799-799`
- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:797-797`
<!-- vulcan:connections:end -->

## Limitations
None beyond being a live environment read intended only for the snapshot builder.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 796.
