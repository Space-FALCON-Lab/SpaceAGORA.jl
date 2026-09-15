---
id: simulation.dynamics_rhs__accumulate_control_effectors_bang
label: _accumulate_control_effectors!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _accumulate_control_effectors!
  lines:
  - 183
  - 183
inputs:
- id: forces
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `forces`.
- id: torques
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `torques`.
- id: rw_torque_body
  type: MVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `rw_torque_body`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: debug_control
  type: Bool
  units: n/a
  required: true
  description: Positional argument `debug_control`.
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
  type: Float64
  units: n/a
  description: Return value of `_accumulate_control_effectors!`; mutates `forces`
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

# _accumulate_control_effectors!

## Purpose
Sums every control effector's force, torque, reaction-wheel torque and mass flow rate for one satellite, the control half of the RHS.

## Design & Implementation
Loops the control effector tuple calling `calcControlForceTorque`, `calcControlMassFlowRate` and `calcReactionWheelTorque`, accumulating into the three `MVector`s and returning the summed mass rate with non-finite contributions dropped. When `debug_control` is set it prints any non-zero force or torque. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `forces` | MVector{3, Float64} | n/a | yes | Positional argument `forces`. |
| in | `torques` | MVector{3, Float64} | n/a | yes | Positional argument `torques`. |
| in | `rw_torque_body` | MVector{3, Float64} | n/a | yes | Positional argument `rw_torque_body`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `debug_control` | Bool | n/a | yes | Positional argument `debug_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_accumulate_control_effectors!`; mutates `forces` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1399-1399`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2116-2116`
- [[simulation.dynamics_rhs_spacecraft_dynamics_fast_control_bang|spacecraft_dynamics_fast_control!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2220-2220`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1747-1747`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:199-199`
- `callees` → [[gnc.momentum_manager_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- `callees` → [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- `callees` → [[gnc.propulsive_maneuvers_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:196-196`
- `callees` → [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:197-197`
- `callees` → [[gnc.robot_arm_control_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- `callees` → [[gnc.robot_arm_control_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:196-196`
- `callees` → [[gnc.rpo_mpc_control_model_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
- `callees` → [[gnc.rpo_mpc_control_model_calccontrolmassflowrate|calcControlMassFlowRate]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:196-196`
- `callees` → [[gnc.rpo_mpc_control_model_calcreactionwheeltorque|calcReactionWheelTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:197-197`
- `callees` → [[gnc.targeting_control_calccontrolforcetorque|calcControlForceTorque]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:195-195`
<!-- vulcan:connections:end -->

## Limitations
The `println` debugging path allocates strings on the hot path when enabled; and dropping non-finite mass rates silently hides an effector bug.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 183.
