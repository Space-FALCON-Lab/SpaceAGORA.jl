---
id: gnc.robot_arm_control_init_robot_arm_joint_mpc
label: init_robot_arm_joint_mpc
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: init_robot_arm_joint_mpc
  lines:
  - 56
  - 56
inputs:
- id: plan
  type: RobotArmPlan
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: dt_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `dt_s` (default `0.1`).
- id: horizon
  type: Integer
  units: n/a
  required: false
  description: Keyword argument `horizon` (default `8`).
- id: q_weight
  type: Real
  units: n/a
  required: false
  description: Keyword argument `q_weight` (default `40.0`).
- id: dq_weight
  type: Real
  units: n/a
  required: false
  description: Keyword argument `dq_weight` (default `6.0`).
- id: torque_weight
  type: Real
  units: n/a
  required: false
  description: Keyword argument `torque_weight` (default `0.1`).
- id: terminal_weight_scale
  type: Real
  units: n/a
  required: false
  description: Keyword argument `terminal_weight_scale` (default `4.0`).
- id: joint_inertia_kg_m2
  type: Any
  units: n/a
  required: false
  description: Keyword argument `joint_inertia_kg_m2` (default `nothing`).
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
  description: Return value of `init_robot_arm_joint_mpc`. Returns `RobotArmJointMPCController(Ad,
    Bd, horizon_i, H, E, F, zeros(nu * horizon_i), in`.
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

# init_robot_arm_joint_mpc

## Purpose
Builds the `RobotArmJointMPCController` for a given plan by discretising each joint as an independent double integrator with its effective inertia and condensing a finite-horizon quadratic tracking cost into the dense matrices `H`, `E` and `F` that the runtime control law solves against.

## Theory & Math
Per-joint discrete dynamics with sample time $\Delta t$ and inertia $J_i$: $\begin{pmatrix} q \\ \dot q\end{pmatrix}_{k+1} = \begin{pmatrix} I & \Delta t I \\ 0 & I\end{pmatrix}\begin{pmatrix} q \\ \dot q\end{pmatrix}_k + \begin{pmatrix} \tfrac{1}{2}\Delta t^2 J^{-1} \\ \Delta t J^{-1}\end{pmatrix}\tau_k$. Stacking $N$ steps, $X = \bar A x_0 + \bar B U$, and the cost $\sum_{k=1}^{N}(x_k - r_k)^T Q_k (x_k - r_k) + \tau_{k-1}^T R \tau_{k-1}$ with $Q_N = s_f Q$ condenses to $\tfrac{1}{2}U^T H U - U^T(F\,\bar r - E x_0)$ where $H = \bar B^T \bar Q \bar B + \bar R$, $E = \bar B^T \bar Q \bar A$, $F = \bar B^T \bar Q$.

## Design & Implementation
Keyword arguments `dt_s=0.1`, `horizon=8`, `q_weight=40.0`, `dq_weight=6.0`, `torque_weight=0.1`, `terminal_weight_scale=4.0` and optional `joint_inertia_kg_m2`. It validates `horizon > 0`, `dt > 0`, inertia length `n = length(plan.q_start)` and positivity, throwing `ArgumentError` otherwise; missing inertia falls back to `_robot_arm_default_joint_inertia(plan)`. The plant is `Ad = [I dt*I; 0 I]`, `Bd = [0.5 dt^2 J^-1; dt J^-1]` with `J^-1 = Diagonal(1 ./ inertia)`. `rpo_prediction_mats(Ad, Bd, horizon)` yields the stacked prediction `Abar`, `Bbar`; weights are `Q = diag(q_weight * 1_n, dq_weight * 1_n)`, `Qf = terminal_weight_scale * Q`, `R = torque_weight * I`, assembled into block-diagonal `Qbar` (with `Qf` in the last block) and `Rbar` by `rpo_block_diag`. Then `H = sym(Bbar' Qbar Bbar + Rbar)`, `E = Bbar' Qbar Abar`, `F = Bbar' Qbar`, and the controller is returned with `U_prev = zeros(n * horizon)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | RobotArmPlan | n/a | yes | Positional argument `plan`. |
| in | `dt_s` | Real | n/a | no | Keyword argument `dt_s` (default `0.1`). |
| in | `horizon` | Integer | n/a | no | Keyword argument `horizon` (default `8`). |
| in | `q_weight` | Real | n/a | no | Keyword argument `q_weight` (default `40.0`). |
| in | `dq_weight` | Real | n/a | no | Keyword argument `dq_weight` (default `6.0`). |
| in | `torque_weight` | Real | n/a | no | Keyword argument `torque_weight` (default `0.1`). |
| in | `terminal_weight_scale` | Real | n/a | no | Keyword argument `terminal_weight_scale` (default `4.0`). |
| in | `joint_inertia_kg_m2` | Any | n/a | no | Keyword argument `joint_inertia_kg_m2` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmJointMPCController | n/a | — | Return value of `init_robot_arm_joint_mpc`. Returns `RobotArmJointMPCController(Ad, Bd, horizon_i, H, E, F, zeros(nu * horizon_i), in`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:69-69`
- `callees` → [[gnc.lqmpc_rpo_block_diag|rpo_block_diag]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:92-92`
- `callees` → [[gnc.lqmpc_rpo_prediction_mats|rpo_prediction_mats]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:88-88`
- `callees` → [[gnc.robot_arm_control__robot_arm_default_joint_inertia|_robot_arm_default_joint_inertia]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:72-72`
- `callees` → [[gnc.robot_arm_control_robotarmjointmpccontroller|RobotArmJointMPCController]] · `callers` · call · `src/gnc/control/robot_arm_control.jl:98-98`
<!-- vulcan:connections:end -->

## Limitations
Joints are modelled as decoupled double integrators, so gravity-free but fully coupled multibody dynamics (Coriolis, inertia coupling between links, base motion) are unmodelled disturbances the MPC must reject through feedback alone. No input or state constraints are included. Because the reference is tracked only at steps 1..N, the current-step reference `ref[:, 1]` is discarded by the control law. The default weights are tuned constants without units documented beyond joint radians and N m.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 56.
