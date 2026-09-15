---
id: simulation.setup__rhs_batch_thread_threshold
label: _rhs_batch_thread_threshold
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_batch_thread_threshold
  lines:
  - 371
  - 371
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
  description: Return value of `_rhs_batch_thread_threshold`.
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

# _rhs_batch_thread_threshold

## Purpose
Minimum number of spacecraft an RHS batch must contain before `:auto` mode will spread the batch across threads, avoiding threading overhead on small constellations.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_BATCH_THREAD_THRESHOLD", 16)`, clamped to at least 1. In `_rhs_batch_parallel_enabled` the check is `num_spacecraft >= env.batch_thread_threshold && Polyester.num_cores() > 1`. Captured once into `RhsPlanEnvConfig.batch_thread_threshold`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_batch_thread_threshold`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:854-854`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:372-372`
<!-- vulcan:connections:end -->

## Limitations
The default 16 is a fixed heuristic that ignores per-satellite cost; a constellation of 8 satellites with expensive harmonics stays serial in `:auto` mode. There is no upper bound. Setting the threshold to 1 with `Polyester.num_cores() > 1` parallelises every batch.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 371.
