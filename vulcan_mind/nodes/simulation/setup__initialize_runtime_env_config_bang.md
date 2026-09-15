---
id: simulation.setup__initialize_runtime_env_config_bang
label: _initialize_runtime_env_config!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_runtime_env_config!
  lines:
  - 912
  - 912
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Nothing
  units: n/a
  description: Return value of `_initialize_runtime_env_config!`; mutates `p` in place.
    Returns `nothing`.
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

# _initialize_runtime_env_config!

## Purpose
Captures the three run-scoped environment snapshots — policy, RHS plan and callback — into shared buffers so hot paths never parse `ENV` again during the run.

## Design & Implementation
Assigns `snapshot_policy_decision_env()`, `_snapshot_rhs_plan_env_config()` and `_snapshot_callback_env_config()` into their `Ref` slots. Called from `run_simulation` setup inside any active engine-config override scope.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_runtime_env_config!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:189-189`

**Downstream**

- `callees` → [[parcore.env_config_snapshot_policy_decision_env|snapshot_policy_decision_env]] · `callers` · call · `src/simulation/engine/setup.jl:913-913`
- `callees` → [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callers` · call · `src/simulation/engine/setup.jl:914-914`
- `callees` → [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callers` · call · `src/simulation/engine/setup.jl:915-915`
<!-- vulcan:connections:end -->

## Limitations
Once captured, changing an environment variable has no effect on the running simulation; hand-built parameters that skip this step fall back to live parsing.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 912.
