---
id: simulation.dynamics_rhs__prepare_rhs_flat_work_packets_bang
label: _prepare_rhs_flat_work_packets!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _prepare_rhs_flat_work_packets!
  lines:
  - 495
  - 495
inputs:
- id: packet_starts
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `packet_starts`.
- id: packet_ends
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `packet_ends`.
- id: packet_costs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `packet_costs`.
- id: packet_elapsed_ns
  type: Vector{Int64}
  units: n/a
  required: true
  description: Positional argument `packet_elapsed_ns`.
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
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
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
  description: Return value of `_prepare_rhs_flat_work_packets!`; mutates `packet_starts`
    in place.
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

# _prepare_rhs_flat_work_packets!

## Purpose
Groups the flat work items into packets of roughly equal estimated cost so the persistent worker pool schedules fewer, larger units and amortises task overhead.

## Design & Implementation
Resizes the four packet vectors if needed, sums estimated item costs, derives a per-packet target from the total, worker count and the minimum packet duration, then walks items in order opening a new packet whenever the running cost would exceed the target. Records each packet's start, end and estimated cost and returns the packet count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `packet_starts` | Vector{Int} | n/a | yes | Positional argument `packet_starts`. |
| in | `packet_ends` | Vector{Int} | n/a | yes | Positional argument `packet_ends`. |
| in | `packet_costs` | Vector{Float64} | n/a | yes | Positional argument `packet_costs`. |
| in | `packet_elapsed_ns` | Vector{Int64} | n/a | yes | Positional argument `packet_elapsed_ns`. |
| in | `work_items` | Vector{Int} | n/a | yes | Positional argument `work_items`. |
| in | `count_items` | Int | n/a | yes | Positional argument `count_items`. |
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_prepare_rhs_flat_work_packets!`; mutates `packet_starts` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1084-1084`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__rhs_flat_item_estimated_cost_ns|_rhs_flat_item_estimated_cost_ns]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:516-516`
- `callees` → [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:518-518`
<!-- vulcan:connections:end -->

## Limitations
Greedy grouping in satellite-major item order, so a single very expensive effector can produce a packet far over target; and elapsed times are recorded per packet, so cost feedback is only as fine as the packets.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 495.
