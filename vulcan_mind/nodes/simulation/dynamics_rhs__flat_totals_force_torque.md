---
id: simulation.dynamics_rhs__flat_totals_force_torque
label: _flat_totals_force_torque
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _flat_totals_force_torque
  lines:
  - 1195
  - 1195
inputs:
- id: totals
  type: Any
  units: n/a
  required: true
  description: Positional argument `totals`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: Any
  units: n/a
  description: Return value of `_flat_totals_force_torque`. Returns `forces, torques`.
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

# _flat_totals_force_torque

## Purpose
Extracts a satellite's summed force and torque from the flat totals matrix into mutable vectors the per-satellite tail can continue accumulating into.

## Design & Implementation
Reads rows one to three into a force `MVector` and rows four to six into a torque `MVector` at column `sat_idx`, returning the pair. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `totals` | Any | n/a | yes | Positional argument `totals`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_flat_totals_force_torque`. Returns `forces, torques`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1350-1350`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The row layout is an implicit contract with every batch kernel and the worker reduction; rows seven and eight of the eight-row totals are reserved and not read here.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1195.
