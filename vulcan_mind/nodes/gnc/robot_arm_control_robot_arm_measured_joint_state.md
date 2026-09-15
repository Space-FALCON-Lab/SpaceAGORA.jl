---
id: gnc.robot_arm_control_robot_arm_measured_joint_state
label: robot_arm_measured_joint_state
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: robot_arm_measured_joint_state
  lines:
  - 147
  - 147
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
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
  description: Return value of `robot_arm_measured_joint_state`. Returns `vcat(q_meas,
    dq_meas)`.
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

# robot_arm_measured_joint_state

## Purpose
Reconstructs the arm's joint angles and joint rates from the simulated link attitudes and angular velocities in the spacecraft state view, providing the feedback measurement `x = [q; dq]` for the joint-space MPC so it tracks the actual arm rather than assuming it is on the reference.

## Design & Implementation
Signature `robot_arm_measured_joint_state(plan::RobotArmPlan, sc_view)`. It returns `nothing` unless `sc_view` has both `arm_q` (4 x n quaternion columns) and `arm_ω` (3 x n body-frame angular velocities). The parent starts as the spacecraft: `parent_q` from `sc_view.q` if present (else `plan.base_pose.quaternion`) and `parent_ω_world = rot(parent_q) * sc_view.ω` if present (else zero). For each joint `i` it forms `child_q`, the relative quaternion `quat_mult(conj(parent_q), child_q)`, the joint axis in parent coordinates `plan.model.joints[i].axis_parent` and in world coordinates `rot(parent_q) * axis_parent`, then sets `q_meas[i] = _robot_arm_control_axis_angle_about(rel_q, axis_parent)` and `dq_meas[i] = dot(child_ω_world - parent_ω_world, axis_world)`. The child becomes the parent for the next joint. Returns `vcat(q_meas, dq_meas)` of length `2n`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `robot_arm_measured_joint_state`. Returns `vcat(q_meas, dq_meas)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:191-191`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- `callees` → [[core.project_unit_quaternion|project_unit_quaternion]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:153-153`
- `callees` → [[core.quaternion_utils_quat_mult|quat_mult]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:159-159`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:155-155`
- `callees` → [[gnc.robot_arm_control__robot_arm_control_axis_angle_about|_robot_arm_control_axis_angle_about]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:163-163`
- `callees` → [[gnc.robot_arm_control__robot_arm_control_quat_conj|_robot_arm_control_quat_conj]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:159-159`
<!-- vulcan:connections:end -->

## Limitations
The chain is assumed serial with joint `i` connecting link `i-1` to link `i`; branched arms are not supported. If `sc_view.q` is absent the base attitude is taken from the plan's nominal pose, silently ignoring any base motion since planning. Measured angles are wrapped to `(-π, π]`, so a joint that has rotated past π reports a discontinuity that the MPC will try to correct with a large torque. The `arm_ω` columns are assumed to be expressed in each link's own body frame.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 147.
