---
id: gnc.robot_arm_control_robotarmjointmpccontroller
label: RobotArmJointMPCController
kind: struct
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: RobotArmJointMPCController
  lines:
  - 9
  - 9
inputs:
- id: Ad
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `Ad`.
- id: Bd
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `Bd`.
- id: horizon
  type: Int
  units: n/a
  required: true
  description: Field `horizon`.
- id: H
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `H`.
- id: E
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `E`.
- id: F
  type: Matrix{Float64}
  units: n/a
  required: true
  description: Field `F`.
- id: U_prev
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `U_prev`.
- id: joint_inertia_kg_m2
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `joint_inertia_kg_m2`.
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
  type: RobotArmJointMPCController
  units: n/a
  description: Constructed `RobotArmJointMPCController`.
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

# RobotArmJointMPCController

## Purpose
Holds the precomputed matrices of a joint-space linear MPC problem for robot-arm trajectory tracking, so that each control update reduces to forming a right-hand side and solving one dense linear system rather than rebuilding the prediction model.

## Design & Implementation
A `mutable struct` with discrete double-integrator matrices `Ad::Matrix{Float64}` (2n x 2n) and `Bd` (2n x n), the prediction `horizon::Int`, the condensed quadratic-cost Hessian `H = Bbar' Qbar Bbar + Rbar` (symmetrised), the state-coupling matrix `E = Bbar' Qbar Abar`, the reference-coupling matrix `F = Bbar' Qbar`, a warm-start buffer `U_prev::Vector{Float64}` of length `n * horizon`, and `joint_inertia_kg_m2::Vector{Float64}`. `init_robot_arm_joint_mpc` builds it; `robot_arm_joint_mpc_control` solves `H U = F r_stack - E x` with `\` and shifts the solution into `U_prev`. `calcControlEffect!` dispatches on `model.controller isa RobotArmJointMPCController`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `Ad` | Matrix{Float64} | n/a | yes | Field `Ad`. |
| in | `Bd` | Matrix{Float64} | n/a | yes | Field `Bd`. |
| in | `horizon` | Int | n/a | yes | Field `horizon`. |
| in | `H` | Matrix{Float64} | n/a | yes | Field `H`. |
| in | `E` | Matrix{Float64} | n/a | yes | Field `E`. |
| in | `F` | Matrix{Float64} | n/a | yes | Field `F`. |
| in | `U_prev` | Vector{Float64} | n/a | yes | Field `U_prev`. |
| in | `joint_inertia_kg_m2` | Vector{Float64} | n/a | yes | Field `joint_inertia_kg_m2`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmJointMPCController | n/a | — | Constructed `RobotArmJointMPCController`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_init_robot_arm_joint_mpc|init_robot_arm_joint_mpc]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:98-98`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The problem is unconstrained: no torque, joint-angle or joint-rate limits enter `H`, `E` or `F`, so the solved torques can exceed actuator capability. `U_prev` is updated but never used as a warm start because the solve is a direct dense factorisation of `H` each call. The dense `H` is `(n*horizon)^2` in size and is refactorised every update; for the default `horizon=8` and a 6-joint arm that is a 48 x 48 solve per step. Storing `Any`-typed `controller` on the effector defeats compile-time specialisation.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 9.
