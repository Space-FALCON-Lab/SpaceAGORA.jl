---
id: simulation.setup__rhs_flat_packet_min_items
label: _rhs_flat_packet_min_items
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_packet_min_items
  lines:
  - 437
  - 437
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
  description: Return value of `_rhs_flat_packet_min_items`.
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

# _rhs_flat_packet_min_items

## Purpose
Minimum number of `(satellite, effector)` work items the flat queue must contain before packet scheduling is considered, since packetising a tiny queue costs more than it saves.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_thread_threshold_env("SPACEAGORA_RHS_FLAT_PACKET_MIN_ITEMS", 128)`, clamped to at least 1. Compared against `num_sats × n_effectors` by the flat-queue planner when the scheduler mode is `:auto`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_rhs_flat_packet_min_items`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:878-878`

**Downstream**

- `callees` → [[parallel.env_config_parse_thread_threshold_env|parse_thread_threshold_env]] · `callers` · call · `src/simulation/engine/setup.jl:438-438`
<!-- vulcan:connections:end -->

## Limitations
The default 128 is tuned for 24+ satellites with 3+ effectors and does not scale with per-item cost; a 64-item queue of expensive harmonics evaluations never gets packet scheduling in `:auto`. No upper bound.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 437.
