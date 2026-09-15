---
id: simulation.setup__effector_allow_with_outer
label: _effector_allow_with_outer
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_allow_with_outer
  lines:
  - 409
  - 409
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
  type: Bool
  units: n/a
  description: Return value of `_effector_allow_with_outer`.
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

# _effector_allow_with_outer

## Purpose
Opt-in that lets the inner effector loop still use threads when an outer parallel layer is active, halving its share budget instead of disabling it, for operators who have verified their thread budget can absorb both.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EFFECTOR_PARALLEL_ALLOW_WITH_OUTER", false)`. Captured into `RhsPlanEnvConfig.effector_allow_with_outer`. In `_dynamic_effector_thread_decision`, `outer_active && allow_with_outer` reduces `share_budget` to `max(1, fld(share_budget, 2))` and the flag is passed through to `thread_policy_decision`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_effector_allow_with_outer`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:858-858`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:410-410`
<!-- vulcan:connections:end -->

## Limitations
Defaults to false because the repository's own thread-allocation handoff notes document livelock-like contention when outer and inner parallelism split the pool; enabling it re-exposes that hazard. The halving is a fixed heuristic, not derived from the actual outer worker count.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 409.
