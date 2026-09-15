---
id: gnc.robot_arm_control_calccontroleffect_bang
label: calcControlEffect!
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: calcControlEffect!
  lines:
  - 183
  - 183
inputs:
- id: model
  type: RobotArmControlEffector
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Nothing
  units: n/a
  description: Return value of `calcControlEffect!`; mutates `model` in place. Returns
    `nothing`.
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

# calcControlEffect!

## Purpose
Controller update hook: at each control tick for the owning spacecraft it samples the arm plan, measures (or assumes) the joint state, runs the joint-space MPC when one is attached, and stores the resulting joint torques in `model.held` for the dynamics to consume until the next update.

## Design & Implementation
Signature `calcControlEffect!(model::RobotArmControlEffector, u, p::ODEParams, t::Float64, sat_idx::Int)`. It returns early unless `sat_idx == model.spacecraft_idx` and `model.plan !== nothing`. Plan time is `plan_t = max(0, t - model.updated_at_s)`; `robot_arm_plan_sample` gives the joint count via `length(sample.q)`, and `joint_torque` starts as zeros. If `model.controller isa RobotArmJointMPCController`, it obtains `sc_view` via `_robot_arm_control_spacecraft_state`, the measurement `x` via `robot_arm_measured_joint_state` (falling back to `_robot_arm_control_reference_state`), the preview `x_ref` via `robot_arm_joint_mpc_reference_preview(plan, plan_t, model.control_dt_s, controller.horizon)`, and writes `robot_arm_joint_mpc_control(controller, x, x_ref)` into `joint_torque`. Finally `model.held` is replaced by a new `RobotArmHeldActuation` carrying `joint_torque_nm` and zero base force and torque. Returns `nothing`; the only mutations are `model.held` and, inside the MPC solve, `controller.U_prev`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RobotArmControlEffector | n/a | yes | Positional argument `model`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `calcControlEffect!`; mutates `model` in place. Returns `nothing`. |
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

- `callees` → [[gnc.robot_arm_control__robot_arm_control_reference_state|_robot_arm_control_reference_state]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:192-192`
- `callees` → [[gnc.robot_arm_control__robot_arm_control_spacecraft_state|_robot_arm_control_spacecraft_state]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:190-190`
- `callees` → [[gnc.robot_arm_control_robot_arm_joint_mpc_reference_preview|robot_arm_joint_mpc_reference_preview]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:193-193`
- `callees` → [[gnc.robot_arm_control_robot_arm_measured_joint_state|robot_arm_measured_joint_state]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:191-191`
- `callees` → [[gnc.robot_arm_control_robotarmheldactuation|RobotArmHeldActuation]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:196-196`
- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:187-187`
- `callees` → [[gncx.robot_arm_control_robot_arm_joint_mpc_control|robot_arm_joint_mpc_control]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:194-194`
<!-- vulcan:connections:end -->

## Limitations
The initial `sample` is computed only to learn the joint count and is otherwise discarded, costing one plan interpolation per tick. Any controller type other than `RobotArmJointMPCController` (including `nothing`) yields silent zero torques. The function does not enforce `control_dt_s` itself; if the caller invokes it more often than that, the MPC is solved at every call with a preview spaced at `control_dt_s`, mismatching the actual update rate. Base force and torque are always zeroed, so the compliance gains on the effector have no path into `held`.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 183.
