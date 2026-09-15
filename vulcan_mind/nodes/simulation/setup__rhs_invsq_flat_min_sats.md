---
id: simulation.setup__rhs_invsq_flat_min_sats
label: _rhs_invsq_flat_min_sats
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_invsq_flat_min_sats
  lines:
  - 962
  - 962
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
  description: Return value of `_rhs_invsq_flat_min_sats`.
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

# _rhs_invsq_flat_min_sats

## Purpose
Reads the minimum satellite count for routing a lone inverse-square or J2 gravity effector through the flat queue, which is cheap enough per satellite that a lower threshold than the general flat path is warranted.

## Design & Implementation
Parses `SPACEAGORA_INVSQ_FLAT_MIN_SATS` through `parse_thread_threshold_env` with a default of 8. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_invsq_flat_min_sats`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1177-1177`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:963-963`
<!-- vulcan:connections:end -->

## Limitations
Unlike the harmonics flat threshold, this one is read live inside `_rhs_execution_plan_uncached` rather than from the `RhsPlanEnvConfig` snapshot, so it costs an environment lookup on every RHS call that reaches that branch.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 962.
