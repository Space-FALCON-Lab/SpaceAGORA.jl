---
id: simulation.dynamics_rhs_constellationinteractionedgeworkitem
label: ConstellationInteractionEdgeWorkItem
kind: struct
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: ConstellationInteractionEdgeWorkItem
  lines:
  - 378
  - 378
inputs:
- id: source_sat_idx
  type: Int
  units: n/a
  required: true
  description: Field `source_sat_idx`.
- id: target_sat_idx
  type: Int
  units: n/a
  required: true
  description: Field `target_sat_idx`.
- id: eff_idx
  type: Int
  units: n/a
  required: true
  description: Field `eff_idx`.
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
  type: ConstellationInteractionEdgeWorkItem
  units: n/a
  description: Constructed `ConstellationInteractionEdgeWorkItem`.
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

# ConstellationInteractionEdgeWorkItem

## Purpose
A work item describing a pairwise interaction between two satellites under one effector, reserved for future inter-satellite effectors such as formation-flying forces or communication constraints.

## Design & Implementation
An immutable struct with `source_sat_idx`, `target_sat_idx` and `eff_idx`. It is the edge counterpart of the integer node work item, which is why `ConstellationExecutionPlan` carries an `edge_count`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source_sat_idx` | Int | n/a | yes | Field `source_sat_idx`. |
| in | `target_sat_idx` | Int | n/a | yes | Field `target_sat_idx`. |
| in | `eff_idx` | Int | n/a | yes | Field `eff_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ConstellationInteractionEdgeWorkItem | n/a | — | Constructed `ConstellationInteractionEdgeWorkItem`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Declared but not yet produced or consumed by any scheduler; it documents intent rather than behaviour.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 378.
