---
id: simulation.setup__policy_env_config
label: _policy_env_config
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _policy_env_config
  lines:
  - 900
  - 900
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_policy_env_config`.
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

# _policy_env_config

## Purpose
Retrieves the run-scoped policy environment snapshot from the parameters, or `nothing` when the run did not capture one.

## Design & Implementation
Resolves shared buffers through `_effector_shared_buffers`, returns the `policy_env_config[]` if the field exists, else `nothing`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_policy_env_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.config__control_callback_thread_decision|_control_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:303-303`
- [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:258-258`
- [[simulation.config__thermal_callback_thread_decision|_thermal_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:334-334`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:633-633`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1051-1051`

**Downstream**

- `callees` → [[simulation.setup__effector_shared_buffers|_effector_shared_buffers]] · `callers` · call · `src/simulation/engine/setup.jl:901-901`
<!-- vulcan:connections:end -->

## Limitations
Callers must handle `nothing` by falling back to live parsing, which the execution planner does.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 900.
