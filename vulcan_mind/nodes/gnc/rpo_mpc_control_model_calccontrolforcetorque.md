---
id: gnc.rpo_mpc_control_model_calccontrolforcetorque
label: calcControlForceTorque
kind: function
source:
  file: src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl
  symbol: calcControlForceTorque
  lines:
  - 39
  - 39
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
  description: Return value of `calcControlForceTorque`. Returns `model.held.force_ii,
    model.held.torque_body`.
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
Effector interface method returning the inertial force and body torque that the rendezvous MPC controller has already computed, so the dynamics right-hand side can apply them at any solver stage.

## Design & Implementation
Dispatches on `RPOMPCControlModel`. If the queried satellite index `i` is not `model.chaser_idx` it returns a pair of zero `SVector{3, Float64}` vectors, ensuring the target and any other spacecraft see no control action. For the chaser it returns `model.held.force_ii` in newtons and `model.held.torque_body` in newton-metres, the values latched by `calcControlEffect!` at the last control tick. This is a zero-order hold: the command is constant between ticks even though the integrator evaluates the derivative many times within a step.

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
| out | `result` | Any | n/a | — | Return value of `calcControlForceTorque`. Returns `model.held.force_ii, model.held.torque_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:362-362`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:386-386`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It is a pure read of mutable model state, so correctness depends entirely on `calcControlEffect!` having run for the current control interval; if the plan buffer was invalid or the controller was `nothing`, `model.held` retains a stale command from an earlier tick and keeps applying it. The held command is not saturated or rate-limited here. Reading shared mutable model state makes the method unsafe if one model instance were evaluated from several threads.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl` line 39.
