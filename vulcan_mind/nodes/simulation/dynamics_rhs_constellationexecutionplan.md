---
id: simulation.dynamics_rhs_constellationexecutionplan
label: ConstellationExecutionPlan
kind: struct
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: ConstellationExecutionPlan
  lines:
  - 368
  - 368
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Field `num_sats`.
- id: n_effectors
  type: Int
  units: n/a
  required: true
  description: Field `n_effectors`.
- id: workers
  type: Int
  units: n/a
  required: true
  description: Field `workers`.
- id: node_count
  type: Int
  units: n/a
  required: true
  description: Field `node_count`.
- id: edge_count
  type: Int
  units: n/a
  required: false
  description: Field `edge_count` (default `0`).
- id: partition
  type: Union{Nothing, Symbol}
  units: n/a
  required: false
  description: Field `partition` (default `nothing`).
- id: use_packets
  type: Bool
  units: n/a
  required: false
  description: Field `use_packets` (default `false`).
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
  description: Constructed `ConstellationExecutionPlan` (keyword constructor via @kwdef).
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

# ConstellationExecutionPlan

## Purpose
Describes one flat-queue RHS evaluation — satellite and effector counts, worker count, work-item count, partition and whether packets are used — so the dispatch and reduction phases agree on the shape of the work.

## Design & Implementation
A `Base.@kwdef struct` with `num_sats`, `n_effectors`, `workers`, `node_count`, `edge_count` defaulting to zero, an optional `partition` symbol and `use_packets` defaulting to false. Built by `_build_constellation_execution_plan!` and adjusted by `_with_packet_scheduler`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Field `num_sats`. |
| in | `n_effectors` | Int | n/a | yes | Field `n_effectors`. |
| in | `workers` | Int | n/a | yes | Field `workers`. |
| in | `node_count` | Int | n/a | yes | Field `node_count`. |
| in | `edge_count` | Int | n/a | no | Field `edge_count` (default `0`). |
| in | `partition` | Union{Nothing, Symbol} | n/a | no | Field `partition` (default `nothing`). |
| in | `use_packets` | Bool | n/a | no | Field `use_packets` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ConstellationExecutionPlan | n/a | — | Constructed `ConstellationExecutionPlan` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__build_constellation_execution_plan_bang|_build_constellation_execution_plan!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:472-472`
- [[simulation.dynamics_rhs__with_packet_scheduler|_with_packet_scheduler]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:484-484`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`edge_count` is reserved for pairwise satellite interaction items that no scheduler yet produces.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 368.
