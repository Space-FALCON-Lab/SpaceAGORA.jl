---
id: simulation.dynamics_rhs__update_rhs_flat_packet_cost_model_bang
label: _update_rhs_flat_packet_cost_model!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _update_rhs_flat_packet_cost_model!
  lines:
  - 617
  - 617
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
- id: packet_count
  type: Int
  units: n/a
  required: true
  description: Positional argument `packet_count`.
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
  type: Nothing
  units: n/a
  description: Return value of `_update_rhs_flat_packet_cost_model!`; mutates `shared_buffers`
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

# _update_rhs_flat_packet_cost_model!

## Purpose
Feeds measured packet execution times back into the per-effector cost estimates so subsequent planning decisions use observed rather than static costs.

## Design & Implementation
For each packet with a positive elapsed time, computes the ratio of measured to estimated packet cost, then for every item in the packet scales that item's effector estimate by the ratio and passes it to `_update_rhs_effector_cost_model!`, which updates the EMA and sample count. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `shared_buffers` | Any | n/a | yes | Positional argument `shared_buffers`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `work_items` | Vector{Int} | n/a | yes | Positional argument `work_items`. |
| in | `packet_starts` | Vector{Int} | n/a | yes | Positional argument `packet_starts`. |
| in | `packet_ends` | Vector{Int} | n/a | yes | Positional argument `packet_ends`. |
| in | `packet_costs` | Vector{Float64} | n/a | yes | Positional argument `packet_costs`. |
| in | `packet_elapsed_ns` | Vector{Int64} | n/a | yes | Positional argument `packet_elapsed_ns`. |
| in | `packet_count` | Int | n/a | yes | Positional argument `packet_count`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_update_rhs_flat_packet_cost_model!`; mutates `shared_buffers` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1178-1178`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:633-633`
- `callees` → [[simulation.dynamics_rhs__rhs_effector_estimated_cost_ns|_rhs_effector_estimated_cost_ns]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:637-637`
- `callees` → [[simulation.dynamics_rhs__rhs_flat_item_eff_idx|_rhs_flat_item_eff_idx]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:636-636`
- `callees` → [[simulation.setup__update_rhs_effector_cost_model_bang|_update_rhs_effector_cost_model!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:638-638`
<!-- vulcan:connections:end -->

## Limitations
Attributes a packet's time proportionally to its items' prior estimates, so within-packet variation between effectors cannot be learned; and the update loops over every item in every packet each RHS call.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 617.
