---
id: gncx.lqmpc_rpo_lqmpc_control
label: rpo_lqmpc_control
kind: function
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: rpo_lqmpc_control
  lines:
  - 140
  - 162
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace supplying the RPO LQ-MPC law with the initialised
    controller, the relative state, and a reference preview over the horizon.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: accel_cmd_rtn
  type: SVector{3,Float64}
  units: m/s^2
  description: First relative acceleration command of the horizon, expressed in the
    radial-transverse-normal frame.
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
# rpo_lqmpc_control

## Purpose
`rpo_lqmpc_control` computes the rendezvous-and-proximity-operations acceleration command by solving a constrained finite-horizon linear quadratic problem and returning only its first move. The state is the six-element relative position and velocity in the radial-transverse-normal frame and the control is the relative acceleration.

## Theory & Math
Over a horizon $N$ the stacked prediction $\bar X = \bar A x_0 + \bar B U$ turns the cost $\sum_k (x_k - r_k)^\top Q (x_k - r_k) + u_k^\top R u_k + (x_N - r_N)^\top Q_f (x_N - r_N)$ into the quadratic program

$$\min_U\; U^\top H U + 2 q^\top U \quad \text{s.t.}\quad u_{\min} \le u_k \le u_{\max},$$

with $H = \tfrac12(\bar B^\top \bar Q \bar B + \bar R)$ symmetrised, $q = E x_0 - F\,\mathrm{vec}(r_{1:N})$, $E = \bar B^\top \bar Q \bar A$ and $F = \bar B^\top \bar Q$. Only $u_0$ is applied.

## Model & Assumptions
The plant is the Hill-Clohessy-Wiltshire linearisation built by `rpo_hcw_continuous_mats` from the target mean motion and discretised by a zero-order hold in `rpo_discretize_zoh`. Prediction matrices are stacked by `rpo_prediction_mats` and the stage and terminal weights are assembled into block-diagonal cost matrices, giving a dense Hessian `H`, a state-coupling matrix `E`, and a reference-coupling matrix `F`. Control bounds enter as a doubled identity inequality matrix `G` with right-hand side `W`, which is why this controller solves a quadratic program instead of a linear system.

## Design & Implementation
`init_rpo_lqmpc` sets up an OSQP model once, with the upper triangle of twice the Hessian as the quadratic term, warm starting enabled, polishing disabled, absolute and relative tolerances of one part in ten thousand, and an iteration cap of one thousand. At each call the linear term is rebuilt as `E * x - F * ref_stack` and pushed with `OSQP.update_q!`, the previous shifted solution is loaded through `OSQP.warm_start_x!`, and the solver runs in place into a preallocated results object. A reference with the wrong row count raises an `ArgumentError`; a reference shorter than the horizon is padded by repeating its last column. Only `:Solved` and `:Solved_inaccurate` statuses are accepted, and the shifted control sequence is written back into `U_prev` for the next warm start.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace supplying the RPO LQ-MPC law with the initialised controller, the relative state, and a reference preview over the horizon. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `accel_cmd_rtn` | SVector{3,Float64} | m/s^2 | — | First relative acceleration command of the horizon, expressed in the radial-transverse-normal frame. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:17-17`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Any other solver status returns a zero acceleration command, so an infeasible or iteration-limited problem is indistinguishable from a genuine coast at the interface. The Hill-Clohessy-Wiltshire model assumes a circular target orbit and small separation, and neglects differential drag and oblateness. Actuator bounds are enforced on the commanded acceleration, not on the thruster forces that ultimately realise it.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl:140-162`, with setup at lines 72-123.
