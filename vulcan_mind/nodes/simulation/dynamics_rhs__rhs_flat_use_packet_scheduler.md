---
id: simulation.dynamics_rhs__rhs_flat_use_packet_scheduler
label: _rhs_flat_use_packet_scheduler
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_flat_use_packet_scheduler
  lines:
  - 565
  - 565
inputs:
- id: shared_buffers
  type: Any
  units: n/a
  required: true
  description: Positional argument `shared_buffers`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: work_items
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `work_items`.
- id: count_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `count_items`.
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
  description: Return value of `_rhs_flat_use_packet_scheduler`.
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

# _rhs_flat_use_packet_scheduler

## Purpose
Decides whether the flat queue should group items into packets, from mode, item count, estimated work and heterogeneity, and the overhead-based self-disable.

## Design & Implementation
Returns false for `:off`, true for `:on` with more than one item, false if self-disabled, and otherwise requires minimum items, work and heterogeneity thresholds.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `work_items` | Vector{Int} | n/a | yes | Positional argument `work_items`. |
| in | `count_items` | Int | n/a | yes | Positional argument `count_items`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_flat_use_packet_scheduler`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1078-1078`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__rhs_flat_packet_work_stats|_rhs_flat_packet_work_stats]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:581-581`
- `callees` → [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:571-571`
<!-- vulcan:connections:end -->

## Limitations
Auto mode depends on the cost model having warmed up.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 565.
