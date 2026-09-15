---
id: gnc.lqmpc_rpo_prediction_mats
label: rpo_prediction_mats
kind: function
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: rpo_prediction_mats
  lines:
  - 35
  - 35
inputs:
- id: Ad
  type: Any
  units: n/a
  required: true
  description: Positional argument `Ad`.
- id: Bd
  type: Any
  units: n/a
  required: true
  description: Positional argument `Bd`.
- id: horizon
  type: Int
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
  description: Return value of `rpo_prediction_mats`. Returns `Abar, Bbar`.
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

# rpo_prediction_mats

## Purpose
Stacks the discrete dynamics over the MPC horizon so the whole predicted state trajectory is one affine function of the initial state and the stacked input sequence.

## Theory & Math
$$
X = \bar{A} x_0 + \bar{B} U,\qquad \bar{A} = \begin{bmatrix} A_d \\ A_d^2 \\ \vdots \\ A_d^N \end{bmatrix},\qquad \bar{B}_{ij} = \begin{cases} A_d^{\,i-j} B_d & j \le i \\ 0 & j > i \end{cases}
$$

where $X$ stacks $x_1 \ldots x_N$, $U$ stacks $u_0 \ldots u_{N-1}$ and $N$ is the horizon.

## Design & Implementation
Precomputes the powers of `Ad` from the identity up to `Ad^horizon` in a vector, then fills `Abar` row-block by row-block with `Ad^i` and `Bbar` as a block lower-triangular matrix whose block at row `i`, column `j` is `Ad^(i-j) * Bd`. Storing the powers once avoids recomputing each product inside the double loop, giving a cost linear in horizon for the powers plus quadratic for the fills.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `Ad` | Any | n/a | yes | Positional argument `Ad`. |
| in | `Bd` | Any | n/a | yes | Positional argument `Bd`. |
| in | `horizon` | Int | n/a | yes | Positional argument `horizon`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_prediction_mats`. Returns `Abar, Bbar`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.lqmpc_init_rpo_lqmpc|init_rpo_lqmpc]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:96-96`
- [[gnc.robot_arm_control_init_robot_arm_joint_mpc|init_robot_arm_joint_mpc]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:88-88`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:42-42`
<!-- vulcan:connections:end -->

## Limitations
Both matrices are dense, so memory grows as `(nx*horizon) * (nu*horizon)`; for long horizons the block structure would be better exploited sparsely, and the caller in `init_rpo_lqmpc` does sparsify the resulting Hessian but not these intermediates.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl` line 35.
