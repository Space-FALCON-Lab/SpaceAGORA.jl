---
id: gnc.robot_arm_control_robot_arm_joint_mpc_reference_preview
label: robot_arm_joint_mpc_reference_preview
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: robot_arm_joint_mpc_reference_preview
  lines:
  - 102
  - 102
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
- id: dt_s
  type: Real
  units: n/a
  required: true
  description: Positional argument `dt_s`.
- id: horizon
  type: Integer
  units: n/a
  required: true
  description: Positional argument `horizon`.
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
  description: Return value of `robot_arm_joint_mpc_reference_preview`. Returns `out`.
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

# robot_arm_joint_mpc_reference_preview

## Purpose
Samples the planned joint trajectory at the MPC's future sample instants to build the reference matrix that `robot_arm_joint_mpc_control` tracks, so the controller anticipates upcoming joint motion rather than reacting only to the current set-point.

## Design & Implementation
Signature `robot_arm_joint_mpc_reference_preview(plan::RobotArmPlan, t_elapsed_s::Real, dt_s::Real, horizon::Integer)`. It allocates `out = zeros(2n, horizon + 1)` with `n = length(plan.q_start)` and, for `j = 0..horizon`, calls `robot_arm_plan_sample(plan, t_elapsed_s + j * dt_s)` and writes `sample.q` into rows `1:n` and `sample.dq` into rows `n+1:2n` of column `j + 1` under `@inbounds`. Column 1 is therefore the current reference and columns 2..horizon+1 are the preview the control law stacks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `t_elapsed_s` | Real | n/a | yes | Positional argument `t_elapsed_s`. |
| in | `dt_s` | Real | n/a | yes | Positional argument `dt_s`. |
| in | `horizon` | Integer | n/a | yes | Positional argument `horizon`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `robot_arm_joint_mpc_reference_preview`. Returns `out`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:193-193`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:106-106`
- `callees` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:106-106`
<!-- vulcan:connections:end -->

## Limitations
Behaviour past the end of the plan depends entirely on `robot_arm_plan_sample`'s extrapolation; this function does not clamp or flag it. The preview is regenerated with a fresh allocation on every control update. Sampling at `j * dt_s` uses the controller's `control_dt_s`, so if the plan was generated at a different resolution the preview interpolates through `robot_arm_plan_sample` rather than reusing plan knots exactly.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 102.
