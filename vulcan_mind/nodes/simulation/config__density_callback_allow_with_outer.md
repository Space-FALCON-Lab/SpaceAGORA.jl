---
id: simulation.config__density_callback_allow_with_outer
label: _density_callback_allow_with_outer
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_callback_allow_with_outer
  lines:
  - 54
  - 54
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
  description: Return value of `_density_callback_allow_with_outer`.
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

# _density_callback_allow_with_outer

## Purpose
Reports whether the density callback may still spawn threads when an outer parallel construct, such as a threaded Monte Carlo or ensemble driver, is already active.

## Design & Implementation
Returns `_parse_bool_env("SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)`. Defaulting to `false` makes nested threading opt-in, since the outer loop normally already saturates the available workers. The snapshot field is passed as `allow_with_outer` to `ParallelPolicy.thread_policy_decision`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_density_callback_allow_with_outer`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.config__thermal_callback_allow_with_outer|_thermal_callback_allow_with_outer]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:128-128`
- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:190-190`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:55-55`
<!-- vulcan:connections:end -->

## Limitations
The flag is a blunt override with no awareness of how many workers the outer construct actually holds, so enabling it can oversubscribe the machine badly. Whether an outer construct is active is itself only a hint obtained from `ParallelPolicy.outer_parallel_active`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 54.
