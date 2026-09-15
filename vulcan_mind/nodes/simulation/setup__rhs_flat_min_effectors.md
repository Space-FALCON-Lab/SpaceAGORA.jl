---
id: simulation.setup__rhs_flat_min_effectors
label: _rhs_flat_min_effectors
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_min_effectors
  lines:
  - 780
  - 780
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
  description: Return value of `_rhs_flat_min_effectors`.
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

# _rhs_flat_min_effectors

## Purpose
Reads the minimum number of dynamic effectors a configuration must carry before the flat effector queue is considered, since with few effectors per-satellite batching already balances well.

## Design & Implementation
Parses `SPACEAGORA_EFFECTOR_FLAT_MIN_EFFECTORS` through `parse_thread_threshold_env` with a default of 3. Declared `@inline`. Captured into the RHS plan snapshot at setup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_flat_min_effectors`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:866-866`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:781-781`
<!-- vulcan:connections:end -->

## Limitations
The planner lets a batch-privileged effector such as harmonics bypass this minimum, so the threshold is not an absolute gate.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 780.
