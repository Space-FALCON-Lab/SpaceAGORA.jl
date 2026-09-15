---
id: simulation.setup__satellite_batch_saturates_pool
label: _satellite_batch_saturates_pool
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _satellite_batch_saturates_pool
  lines:
  - 711
  - 711
inputs:
- id: active_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `active_sats`.
- id: budget
  type: Int
  units: n/a
  required: true
  description: Positional argument `budget`.
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
  description: Return value of `_satellite_batch_saturates_pool`.
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

# _satellite_batch_saturates_pool

## Purpose
Detects the case where satellite-level batching already occupies every thread in the budget, in which case any nested per-effector threading would be pure overhead.

## Design & Implementation
`_satellite_batch_saturates_pool(active_sats::Int, budget::Int)::Bool` returns `active_sats >= budget && budget > 1`. The `budget > 1` clause prevents a single-thread budget from being reported as saturated, which would otherwise mask the serial case. Used by `_rhs_execution_plan_uncached` to force `_with_serial_effector_decision` under `satellite_batch` routing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `active_sats` | Int | n/a | yes | Positional argument `active_sats`. |
| in | `budget` | Int | n/a | yes | Positional argument `budget`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_satellite_batch_saturates_pool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1268-1268`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It counts satellites, not the actual number of Polyester workers the batch loop will spawn, so with chunked scheduling fewer workers may be active than assumed. It also ignores per-satellite cost: many cheap satellites saturate the pool nominally but leave threads idle waiting on the barrier.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 711.
