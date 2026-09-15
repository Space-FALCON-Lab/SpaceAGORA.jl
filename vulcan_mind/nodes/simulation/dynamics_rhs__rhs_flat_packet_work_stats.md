---
id: simulation.dynamics_rhs__rhs_flat_packet_work_stats
label: _rhs_flat_packet_work_stats
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _rhs_flat_packet_work_stats
  lines:
  - 546
  - 546
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_rhs_flat_packet_work_stats`.
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

# _rhs_flat_packet_work_stats

## Purpose
Computes the total estimated work and the cost heterogeneity over the current work items, the two numbers the packet-scheduler decision compares against its thresholds.

## Design & Implementation
Loops the items summing estimated cost and tracking the minimum and maximum; heterogeneity is the maximum over the minimum floored at one nanosecond. Returns the pair.

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
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_rhs_flat_packet_work_stats`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__rhs_flat_use_packet_scheduler|_rhs_flat_use_packet_scheduler]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:581-581`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__rhs_flat_item_estimated_cost_ns|_rhs_flat_item_estimated_cost_ns]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:556-556`
<!-- vulcan:connections:end -->

## Limitations
None beyond inheriting the per-effector granularity of the estimates.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 546.
