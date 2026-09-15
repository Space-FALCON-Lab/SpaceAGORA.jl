---
id: simulation.setup__rhs_plan_step_cache_enabled
label: _rhs_plan_step_cache_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_plan_step_cache_enabled
  lines:
  - 996
  - 996
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
  description: Return value of `_rhs_plan_step_cache_enabled`.
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

# _rhs_plan_step_cache_enabled

## Purpose
Reads whether the per-accepted-step execution-plan cache is on, which lets a multi-stage solver reuse one routing decision across all stages of a step.

## Design & Implementation
Parses `SPACEAGORA_RHS_PLAN_STEP_CACHE` with a default of false. Declared `@inline`. When true, `_rhs_execution_plan` stores its decision in `shared_buffers.rhs_plan_step_cache` and the planet-frame callback clears it once per accepted step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_plan_step_cache_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan|_rhs_execution_plan]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1029-1029`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:997-997`
<!-- vulcan:connections:end -->

## Limitations
It is read on every RHS call by `_rhs_execution_plan`, so the cache's own gate costs an environment lookup per evaluation; extending the snapshot mechanism to this flag would remove that.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 996.
