---
id: simulation.dynamics_rhs__with_packet_scheduler
label: _with_packet_scheduler
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _with_packet_scheduler
  lines:
  - 483
  - 483
inputs:
- id: exec_plan
  type: ConstellationExecutionPlan
  units: n/a
  required: true
  description: Positional argument `exec_plan`.
- id: use_packets
  type: Bool
  units: n/a
  required: true
  description: Positional argument `use_packets`.
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
  description: Return value of `_with_packet_scheduler`.
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

# _with_packet_scheduler

## Purpose
Returns a copy of an execution plan with the `use_packets` flag set to the packet scheduler's decision, since `ConstellationExecutionPlan` is immutable and the decision is made after the plan is built.

## Design & Implementation
Reconstructs the `@kwdef` struct copying `num_sats`, `n_effectors`, `workers`, `node_count`, `edge_count` and `partition` unchanged and substituting `use_packets`. Declared `@inline`, so the copy is a stack-allocated isbits rebuild with no heap traffic.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `exec_plan` | ConstellationExecutionPlan | n/a | yes | Positional argument `exec_plan`. |
| in | `use_packets` | Bool | n/a | yes | Positional argument `use_packets`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ConstellationExecutionPlan | n/a | — | Return value of `_with_packet_scheduler`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1079-1079`

**Downstream**

- `callees` → [[simulation.dynamics_rhs_constellationexecutionplan|ConstellationExecutionPlan]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:484-484`
<!-- vulcan:connections:end -->

## Limitations
Every field must be named explicitly, so adding a field to the plan struct requires updating this copier or the new field silently reverts to its default.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 483.
