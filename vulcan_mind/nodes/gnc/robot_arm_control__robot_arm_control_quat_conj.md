---
id: gnc.robot_arm_control__robot_arm_control_quat_conj
label: _robot_arm_control_quat_conj
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: _robot_arm_control_quat_conj
  lines:
  - 133
  - 133
inputs:
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
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
  type: SVector
  units: n/a
  description: Return value of `_robot_arm_control_quat_conj`. Returns `SVector{4,
    Float64}(-qv[1], -qv[2], -qv[3], qv[4])`.
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

# _robot_arm_control_quat_conj

## Purpose
Returns the conjugate (inverse for unit quaternions) of an attitude quaternion in the scalar-last convention used throughout the robot-arm code, after projecting the input back onto the unit sphere. `robot_arm_measured_joint_state` uses it to express a child link's attitude relative to its parent.

## Design & Implementation
`@inline function _robot_arm_control_quat_conj(q)`. It calls `project_unit_quaternion(q)` to renormalise (guarding against integrator drift) and returns `SVector{4,Float64}(-qv[1], -qv[2], -qv[3], qv[4])`, negating the vector part and keeping the scalar part `qv[4]`. Combined with `quat_mult`, `quat_mult(conj(parent), child)` yields the relative rotation from parent to child frame.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `_robot_arm_control_quat_conj`. Returns `SVector{4, Float64}(-qv[1], -qv[2], -qv[3], qv[4])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_robot_arm_measured_joint_state|robot_arm_measured_joint_state]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:159-159`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- `callees` → [[core.project_unit_quaternion|project_unit_quaternion]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:134-134`
<!-- vulcan:connections:end -->

## Limitations
The scalar-last layout is assumed implicitly; passing a scalar-first quaternion produces a wrong but silently valid result. Renormalisation hides, rather than reports, a badly non-unit input such as an all-zero quaternion, whose projection behaviour is defined by `project_unit_quaternion`. The function accepts any 4-element indexable `q` without a length check.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 133.
