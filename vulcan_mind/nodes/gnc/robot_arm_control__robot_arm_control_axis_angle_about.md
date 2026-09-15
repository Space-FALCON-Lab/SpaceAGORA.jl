---
id: gnc.robot_arm_control__robot_arm_control_axis_angle_about
label: _robot_arm_control_axis_angle_about
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: _robot_arm_control_axis_angle_about
  lines:
  - 139
  - 139
inputs:
- id: q_rel
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_rel`.
- id: axis
  type: Any
  units: n/a
  required: true
  description: Positional argument `axis`.
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
  description: Return value of `_robot_arm_control_axis_angle_about`. Returns `2.0
    * atan(dot(SVector{3, Float64}(q[1], q[2], q[3]), a), q[4])`.
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

# _robot_arm_control_axis_angle_about

## Purpose
Extracts the signed rotation angle of a relative quaternion about a specified joint axis, giving the measured joint angle for a revolute joint from parent and child link attitudes.

## Theory & Math
For a unit quaternion $q = (\mathbf{v}, w)$ (scalar last) and unit axis $\mathbf{a}$, the twist angle about $\mathbf{a}$ is $\theta = 2\,\operatorname{atan2}(\mathbf{v}\cdot\mathbf{a},\ w)$, which equals the full rotation angle when the rotation is purely about $\mathbf{a}$ and otherwise returns the component of rotation about that axis.

## Design & Implementation
`@inline function _robot_arm_control_axis_angle_about(q_rel, axis)`. It projects `q_rel` to a unit quaternion, flips its sign when the scalar part `q[4] < 0` so the angle lies in the short-arc branch, converts `axis` to an `SVector{3,Float64}`, and returns `2 * atan(dot(q_vec, a), q[4])` where `q_vec = (q[1], q[2], q[3])`. This is the swing-twist decomposition's twist angle about `a`, using the two-argument `atan` for a full-range signed result.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_rel` | Any | n/a | yes | Positional argument `q_rel`. |
| in | `axis` | Any | n/a | yes | Positional argument `axis`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_robot_arm_control_axis_angle_about`. Returns `2.0 * atan(dot(SVector{3, Float64}(q[1], q[2], q[3]), a), q[4])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_robot_arm_measured_joint_state|robot_arm_measured_joint_state]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:163-163`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- `callees` → [[core.project_unit_quaternion|project_unit_quaternion]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:140-140`
<!-- vulcan:connections:end -->

## Limitations
`axis` is not normalised, so a non-unit joint axis scales the projected sine term and biases the angle. Any swing component (rotation about axes perpendicular to `a`, which a rigid revolute joint should not exhibit but numerical drift can introduce) is folded into the result without diagnostic. The sign flip on `q[4]` limits the recoverable angle to `(-π, π]`, so multi-turn joint angles wrap.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 139.
