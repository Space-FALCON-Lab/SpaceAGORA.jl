---
id: gnc.robot_arm_control__robot_arm_control_reference_state
label: _robot_arm_control_reference_state
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: _robot_arm_control_reference_state
  lines:
  - 172
  - 172
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: t_elapsed_s
  type: Real
  units: n/a
  required: true
  description: Positional argument `t_elapsed_s`.
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
  description: Return value of `_robot_arm_control_reference_state`. Returns `vcat(sample.q,
    sample.dq)`.
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

# _robot_arm_control_reference_state

## Purpose
Returns the planned joint state vector `[q; dq]` at an elapsed plan time, used as a stand-in for the measured state when the spacecraft view lacks arm attitude and rate fields (for example in reduced-state simulations or before the arm states are allocated).

## Design & Implementation
Signature `_robot_arm_control_reference_state(plan::RobotArmPlan, t_elapsed_s::Real)`. It calls `robot_arm_plan_sample(plan, t_elapsed_s)` and concatenates `sample.q` and `sample.dq` with `vcat`, giving a `Vector{Float64}` of length `2n`. `calcControlEffect!` uses it as the fallback for `x` when `robot_arm_measured_joint_state` returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `t_elapsed_s` | Real | n/a | yes | Positional argument `t_elapsed_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_control_reference_state`. Returns `vcat(sample.q, sample.dq)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:192-192`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:173-173`
<!-- vulcan:connections:end -->

## Limitations
Feeding the reference back as the measurement turns the MPC into open-loop feedforward: the tracking error term is identically zero, so any real disturbance on the arm goes uncorrected and the caller receives no warning that feedback is absent. Behaviour for negative or past-end `t_elapsed_s` is whatever `robot_arm_plan_sample` implements. A fresh vector is allocated on each call.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 172.
