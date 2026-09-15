---
id: simulation.setup__rhs_batch_parallel_mode
label: _rhs_batch_parallel_mode
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_batch_parallel_mode
  lines:
  - 367
  - 367
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
  description: Return value of `_rhs_batch_parallel_mode`.
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

# _rhs_batch_parallel_mode

## Purpose
Selects whether the outer satellite batch of an RHS evaluation is distributed across threads, independently of the inner per-effector parallelism.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_RHS_BATCH_PARALLEL"; default="auto")` as `:off`, `:on`, or `:auto`. Stored in `RhsPlanEnvConfig.batch_parallel_mode` by `_snapshot_rhs_plan_env_config` and interpreted by `_rhs_batch_parallel_enabled`, where `:auto` defers to the satellite-count threshold and core count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_rhs_batch_parallel_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:853-853`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/engine/setup.jl:368-368`
<!-- vulcan:connections:end -->

## Limitations
Even `:on` is overridden when `profile_forces_serial` is set, so the mode is not the final word. Malformed values throw at snapshot time, which is during `run_simulation` setup rather than at process start.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 367.
