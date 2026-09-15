---
id: simulation.setup__rhs_execution_mode_env
label: _rhs_execution_mode_env
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_execution_mode_env
  lines:
  - 758
  - 758
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
  description: Return value of `_rhs_execution_mode_env`.
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

# _rhs_execution_mode_env

## Purpose
Parses the operator's requested RHS execution strategy (automatic routing, serial, satellite batch, per-satellite effector reduce, or flat constellation effector queue) into the symbol the planner dispatches on.

## Design & Implementation
Reads `SPACEAGORA_RHS_EXECUTION_MODE` via `_engine_env_get` with default `"auto"`, lowercases and strips it. `"auto"` or empty maps to `:auto`; `"serial"`, `"off"`, `"none"` to `:serial`; `"satellite"`, `"satellite_batch"`, `"batch"` to `:satellite_batch`; `"per_satellite"`, `"per_satellite_effector_reduce"`, `"effector"` to `:per_satellite_effector_reduce`; `"flat"`, `"flat_constellation"`, `"flat_constellation_effector_queue"` to `:flat_constellation_effector_queue`. Anything else throws `ArgumentError` listing the five canonical names. It is the first field captured by `_snapshot_rhs_plan_env_config`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_rhs_execution_mode_env`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:851-851`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/setup.jl:759-759`
<!-- vulcan:connections:end -->

## Limitations
Forcing a non-auto mode bypasses the eligibility checks that protect against unsupported effector combinations, so `flat` with a non-thread-safe effector is rejected later (or not at all) rather than here. The synonym table is long and only partially mirrored in the error message.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 758.
