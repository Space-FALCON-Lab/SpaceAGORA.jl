---
id: gnc.targeting_control_calccontrolforcetorque
label: calcControlForceTorque
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: calcControlForceTorque
  lines:
  - 272
  - 272
inputs:
- id: model
  type: AerobrakingEnergyDepletionControlModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `calcControlForceTorque`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# calcControlForceTorque

## Purpose
Satisfies the control-effector interface for both the energy-depletion controller and the panel effector by reporting no direct force or torque.

## Design & Implementation
Two methods, one per model type, each returning a pair of zero `SVector{3}`. The controllers act by changing panel geometry, and the resulting aerodynamic force is produced by the drag effector, so contributing a direct wrench here would double count.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerobrakingEnergyDepletionControlModel | n/a | yes | Positional argument `model`. |
| in | `u` | AbstractVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `calcControlForceTorque`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:386-386`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The interface has no way to signal that the effector acts through geometry, so an engine diagnostic listing control effectors by their force contribution reports these as inert.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 272.
