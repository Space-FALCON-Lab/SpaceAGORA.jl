---
id: gncx.robot_arm_control_robot_arm_joint_mpc_control
label: robot_arm_joint_mpc_control
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: robot_arm_joint_mpc_control
  lines:
  - 114
  - 130
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the robot-arm joint MPC law with a
    preconstructed controller, the measured joint state, and a joint reference preview.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: joint_accel_cmd
  type: Vector{Float64}
  units: rad/s^2
  description: First joint acceleration command of the finite-horizon solution, applied
    until the next control update.
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
# robot_arm_joint_mpc_control

## Purpose
`robot_arm_joint_mpc_control` evaluates the robot-arm joint model predictive control law and returns only the first control move of the horizon, in the usual receding-horizon fashion. Its state is the stacked joint positions and joint rates and its control is joint acceleration.

## Model & Assumptions
The controller is an unconstrained finite-horizon linear quadratic problem whose prediction, cost, and gain matrices are precomputed once by `init_robot_arm_joint_mpc` into the `RobotArmJointMPCController` fields `Ad`, `Bd`, `H`, `E`, and `F`. Because there are no inequality constraints, the optimum is the solution of a single dense linear system rather than a quadratic program, which is what distinguishes this law from the RPO controller in the same package. The reference preview is produced separately by `robot_arm_joint_mpc_reference_preview`, which samples the planned arm trajectory at the control step across the horizon and stacks positions above rates.

## Design & Implementation
The method validates that the reference has exactly as many rows as the state dimension and throws an `ArgumentError` naming both the received and expected row counts when it does not. A reference shorter than the horizon is not rejected: it is padded by repeating its final column, so a plan that ends before the horizon is treated as a hold at the terminal pose. The reference columns from the second onward are flattened into `ref_stack`, the right-hand side is formed as `F * ref_stack - E * x`, and the solve `H \ rhs` yields the whole control sequence. The sequence is then shifted into `U_prev` with the last move duplicated, so the next call starts from a warm sequence even though this law does not need one to converge.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the robot-arm joint MPC law with a preconstructed controller, the measured joint state, and a joint reference preview. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `joint_accel_cmd` | Vector{Float64} | rad/s^2 | — | First joint acceleration command of the finite-horizon solution, applied until the next control update. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:194-194`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The law is unconstrained, so joint torque, rate, and travel limits are not enforced and must be respected by the plan and the inertia model instead. The plant is the double-integrator joint model built at initialisation and ignores arm-base coupling, gravity gradient, and joint flexibility. A singular or ill-conditioned `H` propagates as a solve failure rather than a graceful degradation.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl:114-130`.
