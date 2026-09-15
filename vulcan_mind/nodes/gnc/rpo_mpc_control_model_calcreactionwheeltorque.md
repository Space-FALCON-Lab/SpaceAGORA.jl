---
id: gnc.rpo_mpc_control_model_calcreactionwheeltorque
label: calcReactionWheelTorque
kind: function
source:
  file: src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl
  symbol: calcReactionWheelTorque
  lines:
  - 45
  - 45
inputs:
- id: model
  type: RPOMPCControlModel
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
  type: Any
  units: n/a
  description: Return value of `calcReactionWheelTorque`. Returns `model.held.rw_torque_body`.
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

# calcReactionWheelTorque

## Purpose
Reports the portion of the held body torque attributable to reaction wheels, so that wheel momentum accounting stays separate from thruster-generated torque.

## Design & Implementation
Returns `nothing` when the queried index `i` differs from `model.chaser_idx`, signalling that the spacecraft has no wheel contribution, and otherwise returns `model.held.rw_torque_body` in newton-metres. Like `calcControlForceTorque`, the value is held constant between control ticks, since `calcControlEffect!` writes it once per tick from `rpo_reaction_wheel_torque_command`. Note that the same wheel torque is also summed into `model.held.torque_body`, so the total body torque already includes it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RPOMPCControlModel | n/a | yes | Positional argument `model`. |
| in | `u` | AbstractVector | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `calcReactionWheelTorque`. Returns `model.held.rw_torque_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolmassflowrate|calcControlMassFlowRate]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:345-345`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:197-197`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Returning `nothing` rather than a zero vector forces every caller to branch, and a caller that forwards the result into arithmetic without a null check raises. There is no wheel speed, saturation or momentum-dump logic at this level, so a commanded torque that the wheels could not physically deliver is still reported verbatim. The value is stale whenever the control tick was skipped.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl` line 45.
