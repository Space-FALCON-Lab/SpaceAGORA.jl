---
id: gncx.rpo_mpc_control_model_calccontroleffect_bang
label: calcControlEffect!
kind: function
source:
  file: src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl
  symbol: calcControlEffect!
  lines:
  - 7
  - 36
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace dispatching the RPO control update for the chaser
    spacecraft at the current simulation time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: held_actuation
  type: RPOHeldActuation
  units: mixed
  description: Cached inertial force, body torque, thruster forces, and wheel torque
    held until the next control update.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# calcControlEffect!

## Purpose
This `calcControlEffect!` method is the RPO control update: it runs once per control tick for the chaser spacecraft, closes the loop from relative state to actuator commands, and caches the result so the dynamics can read it at every integrator stage.

## Model & Assumptions
Three guards precede any computation. The call must be for the chaser index, the plan buffer must be marked valid, and the controller must have been initialised; any failure returns without touching the held actuation, leaving the previous command in force. The relative state is formed by `inertial_to_rtn_relative_state` from the chaser and target position and velocity pairs, converting the pair of inertial states into the radial-transverse-normal frame in which the Hill-Clohessy-Wiltshire model is written. Plan time is measured as the elapsed interval since the plan buffer was published, floored at zero, which keeps the reference indexing correct when a plan is republished mid-passage.

## Design & Implementation
The chain from command to wrench is explicit. `rpo_ref_preview` samples the plan over the controller horizon, `rpo_lqmpc_control` returns a relative acceleration in the rotating frame, `rtn_accel_to_inertial` maps it back to inertial axes, and multiplication by the instantaneous chaser mass converts acceleration to force. That inertial force is rotated into the body frame through the attitude quaternion, allocated across the six thrusters by `rpo_allocate_six_axis_thrusters`, and the realised wrench is recovered by `rpo_thruster_wrench_body`. Because allocation saturates against per-thruster limits, the realised force is generally not the desired force, and the routine deliberately rotates the realised body force back to inertial rather than reusing the desired one. Thruster and wheel torques are summed into a single body torque, and the four quantities are stored in a freshly constructed `RPOHeldActuation`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace dispatching the RPO control update for the chaser spacecraft at the current simulation time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `held_actuation` | RPOHeldActuation | mixed | — | Cached inertial force, body torque, thruster forces, and wheel torque held until the next control update. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:377-377`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:19-19`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:23-23`
- `callees` → [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:14-14`
- `callees` → [[core.reference_system_rtn_accel_to_inertial|rtn_accel_to_inertial]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:18-18`
- `callees` → [[gnc.lqmpc_rpo_ref_preview|rpo_ref_preview]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:16-16`
- `callees` → [[gnc.rpo_control_types_rpoheldactuation|RPOHeldActuation]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:29-29`
- `callees` → [[gnc.rpo_mpc_control_model__rpo_control_state_pos_vel|_rpo_control_state_pos_vel]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:12-12`
- `callees` → [[gnc.thruster_allocator_rpo_thruster_wrench_body|rpo_thruster_wrench_body]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:26-26`
- `callees` → [[gncx.lqmpc_rpo_lqmpc_control|rpo_lqmpc_control]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:17-17`
- `callees` → [[gncx.reaction_wheel_allocator_rpo_reaction_wheel_torque_command|rpo_reaction_wheel_torque_command]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:27-27`
- `callees` → [[gncx.thruster_allocator_rpo_allocate_six_axis_thrusters|rpo_allocate_six_axis_thrusters]] · `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:25-25`
<!-- vulcan:connections:end -->

## Limitations
The realised force is fed forward without any anti-windup or reference adjustment, so persistent thruster saturation degrades tracking without informing the controller. Mass is read from the state at the control tick and held constant until the next one. No plume impingement, thruster minimum impulse bit, or actuator lag is modelled, and the method silently does nothing for any spacecraft other than the configured chaser.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:7-36`.
