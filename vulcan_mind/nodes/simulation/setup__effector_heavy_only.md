---
id: simulation.setup__effector_heavy_only
label: _effector_heavy_only
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_heavy_only
  lines:
  - 413
  - 413
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
  description: Return value of `_effector_heavy_only`.
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

# _effector_heavy_only

## Purpose
When true, `:auto` mode threads the effector loop only if the estimated work per worker exceeds the nanosecond target, preventing threading of cheap force evaluations where task overhead dominates.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY", true)`; on by default. Stored as `RhsPlanEnvConfig.effector_heavy_only` and forwarded as the `heavy_only` keyword to `thread_policy_decision`, alongside the computed `heavy_work` boolean.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_effector_heavy_only`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:859-859`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:414-414`
<!-- vulcan:connections:end -->

## Limitations
Disabling it lets item count alone (via `_effector_thread_threshold`, default 2) trigger threading, which the source comments note measured 1.5 to 2× slower for small workloads. The flag does not affect `:on` mode.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 413.
