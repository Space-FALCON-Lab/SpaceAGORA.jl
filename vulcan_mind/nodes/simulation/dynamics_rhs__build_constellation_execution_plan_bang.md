---
id: simulation.dynamics_rhs__build_constellation_execution_plan_bang
label: _build_constellation_execution_plan!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _build_constellation_execution_plan!
  lines:
  - 456
  - 456
inputs:
- id: work_items
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `work_items`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
- id: partition
  type: Union{Nothing, Symbol}
  units: n/a
  required: true
  description: Positional argument `partition`.
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
  type: ConstellationExecutionPlan
  units: n/a
  description: Return value of `_build_constellation_execution_plan!`; mutates `work_items`
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

# _build_constellation_execution_plan!

## Purpose
Assembles the `ConstellationExecutionPlan` for one flat-queue RHS call after preparing the work-item list, giving the scheduler a single record of what it must dispatch.

## Design & Implementation
Calls `_prepare_rhs_flat_work_items!` to fill `work_items` and obtain the node count, then constructs the plan with satellite count, effector count, workers, node count, zero interaction edges, the partition and `use_packets=false`. The packet decision is made afterwards by `_with_packet_scheduler`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `work_items` | Vector{Int} | n/a | yes | Positional argument `work_items`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `partition` | Union{Nothing, Symbol} | n/a | yes | Positional argument `partition`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ConstellationExecutionPlan | n/a | — | Return value of `_build_constellation_execution_plan!`; mutates `work_items` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1067-1067`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__prepare_rhs_flat_work_items_bang|_prepare_rhs_flat_work_items!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:465-465`
- `callees` → [[simulation.dynamics_rhs_constellationexecutionplan|ConstellationExecutionPlan]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:472-472`
<!-- vulcan:connections:end -->

## Limitations
`edge_count` is always zero because pairwise satellite interaction items are declared but no scheduler produces them yet.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 456.
