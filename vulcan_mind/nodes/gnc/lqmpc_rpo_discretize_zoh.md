---
id: gnc.lqmpc_rpo_discretize_zoh
label: rpo_discretize_zoh
kind: function
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: rpo_discretize_zoh
  lines:
  - 24
  - 24
inputs:
- id: A
  type: Any
  units: n/a
  required: true
  description: Positional argument `A`.
- id: B
  type: Any
  units: n/a
  required: true
  description: Positional argument `B`.
- id: dt
  type: Real
  units: n/a
  required: true
  description: Positional argument `dt`.
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
  description: Return value of `rpo_discretize_zoh`. Returns `Md[1:nx, 1:nx], Md[1:nx,
    nx+1:end]`.
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

# rpo_discretize_zoh

## Purpose
Converts continuous linear dynamics into the discrete-time pair the MPC predicts over, assuming the control input is held constant across each step.

## Theory & Math
$$
\exp\left(\begin{bmatrix} A & B \\ 0 & 0 \end{bmatrix} \Delta t\right) = \begin{bmatrix} A_d & B_d \\ 0 & I \end{bmatrix}
$$

so that $x_{k+1} = A_d x_k + B_d u_k$ exactly when $u$ is constant over each interval of length $\Delta t$.

## Design & Implementation
Forms the augmented square matrix of size `nx + nu` with `A` in the top-left block and `B` in the top-right, zero elsewhere, and takes its matrix exponential scaled by `dt` through `exp`. The top-left block of the result is the discrete state matrix and the top-right block the discrete input matrix. This single-exponential construction is exact for zero-order hold and avoids evaluating the input integral separately.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `A` | Any | n/a | yes | Positional argument `A`. |
| in | `B` | Any | n/a | yes | Positional argument `B`. |
| in | `dt` | Real | n/a | yes | Positional argument `dt`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_discretize_zoh`. Returns `Md[1:nx, 1:nx], Md[1:nx, nx+1:end]`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.lqmpc_init_rpo_lqmpc|init_rpo_lqmpc]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:91-91`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:30-30`
<!-- vulcan:connections:end -->

## Limitations
A dense matrix exponential is computed on every call, which is acceptable at initialisation but not per control tick; the caller is expected to cache the result. No check is made that `dt` is positive.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl` line 24.
