---
id: simulation.setup__rhs_flat_min_sats
label: _rhs_flat_min_sats
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_min_sats
  lines:
  - 776
  - 776
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
  description: Return value of `_rhs_flat_min_sats`.
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

# _rhs_flat_min_sats

## Purpose
Reads the minimum active satellite count before the flat effector queue is considered, reflecting that flat scheduling only pays off once the satellite-by-effector work grid is large.

## Design & Implementation
Parses `SPACEAGORA_EFFECTOR_FLAT_MIN_SATS` through `parse_thread_threshold_env` with a default of 24. Declared `@inline`. Captured into `RhsPlanEnvConfig.flat_min_sats` at setup and consulted by both the harmonics flat route and the general flat route.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_flat_min_sats`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:865-865`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:777-777`
<!-- vulcan:connections:end -->

## Limitations
The default was tuned on constellation-scaling benchmarks and is high for a small multi-satellite study, which will therefore never see the flat path without an override.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 776.
