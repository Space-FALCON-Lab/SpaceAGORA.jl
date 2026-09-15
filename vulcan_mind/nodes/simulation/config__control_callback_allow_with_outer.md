---
id: simulation.config__control_callback_allow_with_outer
label: _control_callback_allow_with_outer
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _control_callback_allow_with_outer
  lines:
  - 106
  - 106
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
  description: Return value of `_control_callback_allow_with_outer`.
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

# _control_callback_allow_with_outer

## Purpose
Reports whether the control callback may thread while an outer parallel construct is already running.

## Design & Implementation
Returns `_parse_bool_env("SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER", false)`, so nested threading is opt-in. The snapshot field `control_allow_with_outer` reaches `ParallelPolicy.thread_policy_decision` as its `allow_with_outer` keyword, where it is weighed against the `outer_active` hint.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_control_callback_allow_with_outer`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:199-199`

**Downstream**

- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:107-107`
<!-- vulcan:connections:end -->

## Limitations
Enabling it can oversubscribe the process because neither this flag nor the policy layer knows how many workers the outer construct holds. The outer-activity signal itself is only a hint, so a false negative there lets nesting happen even with the flag left at its default.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 106.
