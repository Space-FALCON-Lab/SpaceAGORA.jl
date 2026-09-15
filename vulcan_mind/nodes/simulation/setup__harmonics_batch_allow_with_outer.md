---
id: simulation.setup__harmonics_batch_allow_with_outer
label: _harmonics_batch_allow_with_outer
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _harmonics_batch_allow_with_outer
  lines:
  - 836
  - 836
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
  description: Return value of `_harmonics_batch_allow_with_outer`.
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

# _harmonics_batch_allow_with_outer

## Purpose
Reads whether the harmonics batch and flat-queue RHS modes may run while an outer campaign-level parallel region is already active.

## Design & Implementation
Parses `SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER` through `parse_bool_env` with a default of false. `@inline`. The default of false prevents inner and outer threading from oversubscribing the pool.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_harmonics_batch_allow_with_outer`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:874-874`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:837-837`
- `callees` → [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callers` · feedback · `src/simulation/engine/setup.jl:841-841`
<!-- vulcan:connections:end -->

## Limitations
Live environment read; the snapshotted `RhsPlanEnvConfig.harmonics_batch_allow_with_outer` is what the hot path actually consults.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 836.
