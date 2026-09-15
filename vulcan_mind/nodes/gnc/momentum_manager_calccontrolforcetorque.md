---
id: gnc.momentum_manager_calccontrolforcetorque
label: calcControlForceTorque
kind: function
source:
  file: src/gnc/control/momentum_manager.jl
  symbol: calcControlForceTorque
  lines:
  - 104
  - 104
inputs:
- id: model
  type: MagneticMomentumManagerModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: Any
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
  type: Any
  units: n/a
  description: Return value of `calcControlForceTorque`. Returns `zero3, model.held_torque_body`.
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
Presents the momentum manager's output through the standard control-effector interface, contributing torque but never force.

## Design & Implementation
Returns a pair of zero vectors when the satellite index does not match the model's own. For the matching satellite it returns zero force and `held_torque_body`, the torque latched by the most recent `calcControlEffect!` call. Force is structurally zero because magnetorquers exert a couple on the vehicle rather than a net translational thrust.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | MagneticMomentumManagerModel | n/a | yes | Positional argument `model`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `calcControlForceTorque`. Returns `zero3, model.held_torque_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/momentum_manager.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:386-386`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the torque is whatever was latched previously, calling this without an intervening effect update returns a stale command; the interface offers no way to signal that staleness to the caller.

## Provenance
Mapped from `src/gnc/control/momentum_manager.jl` line 104.
