---
id: dynamics.aerodynamic_wrench_models__multibody_parallel_mode
label: _multibody_parallel_mode
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _multibody_parallel_mode
  lines:
  - 8
  - 8
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
  description: Return value of `_multibody_parallel_mode`.
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

# _multibody_parallel_mode

## Purpose
Reads the process-wide multibody threading mode symbol from the `SPACEAGORA_MULTIBODY_PARALLEL` environment variable through the shared `ParallelPolicy` parser.

## Design & Implementation
A one-line `@inline` wrapper returning `ParallelPolicy.parse_parallel_mode_env("SPACEAGORA_MULTIBODY_PARALLEL")::Symbol`. It is consulted by `_multibody_thread_decision` on every aerodynamic force evaluation and the resulting mode is forwarded to `ParallelPolicy.thread_policy_decision` and recorded in policy observations.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_multibody_parallel_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:36-36`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[parallel.env_config_parse_parallel_mode_env|parse_parallel_mode_env]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:9-9`
<!-- vulcan:connections:end -->

## Limitations
The environment variable is parsed on every call rather than cached, so the string lookup and parse sit inside the RHS hot path. The accepted vocabulary and default are defined entirely in `ParallelPolicy`, so this file cannot validate the value on its own.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 8.
