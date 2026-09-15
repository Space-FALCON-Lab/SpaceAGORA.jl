---
id: simulation.setup__effector_parallel_mode
label: _effector_parallel_mode
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_parallel_mode
  lines:
  - 363
  - 363
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
  description: Return value of `_effector_parallel_mode`.
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

# _effector_parallel_mode

## Purpose
Reads the operator's choice of whether per-effector force evaluation inside a single RHS call may run across threads: forced off, forced on, or decided automatically from cost estimates.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_EFFECTOR_PARALLEL"; default="auto")`, yielding `:off`, `:on`, or `:auto`. Synonyms such as `serial`/`threads` are accepted; anything else throws `ArgumentError`. The value is snapshotted into the RHS environment configuration and consumed by `_dynamic_effector_thread_decision`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_effector_parallel_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:855-855`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/simulation/engine/setup.jl:364-364`
<!-- vulcan:connections:end -->

## Limitations
`:on` does not bypass `_dynamic_effectors_parallel_supported`; effectors flagged as not thread-safe still run serially, so forcing on can be silently ineffective. The variable is read from raw `ENV`, not the engine override layer.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 363.
