---
id: dynamics.aerodynamic_wrench_models__multibody_thread_threshold
label: _multibody_thread_threshold
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _multibody_thread_threshold
  lines:
  - 12
  - 12
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
  description: Return value of `_multibody_thread_threshold`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _multibody_thread_threshold

## Purpose
Returns the minimum number of spacecraft links required before the aerodynamic loop is allowed to use threads, read from `SPACEAGORA_MULTIBODY_THREAD_THRESHOLD` with default 4.

## Design & Implementation
Delegates to `ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD", 4)` and returns an `Int`. The value is passed as `threshold` into `ParallelPolicy.thread_policy_decision` by `_multibody_thread_decision`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_multibody_thread_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:37-37`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
Default `4` is a hard-coded literal duplicated with `_multibody_max_threads`. Non-integer or negative environment values are handled by the parser's rules, not here. Re-parsed on every RHS call.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 12.
