---
id: gnc.lqmpc_init_rpo_lqmpc
label: init_rpo_lqmpc
kind: function
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: init_rpo_lqmpc
  lines:
  - 89
  - 89
inputs:
- id: n
  type: Any
  units: n/a
  required: true
  description: Positional argument `n`.
- id: dt
  type: Any
  units: n/a
  required: true
  description: Positional argument `dt`.
- id: Q
  type: Any
  units: n/a
  required: true
  description: Positional argument `Q`.
- id: R
  type: Any
  units: n/a
  required: true
  description: Positional argument `R`.
- id: Qf
  type: Any
  units: n/a
  required: true
  description: Positional argument `Qf`.
- id: horizon
  type: Any
  units: n/a
  required: true
  description: Positional argument `horizon`.
- id: u_min
  type: Any
  units: n/a
  required: false
  description: Keyword argument `u_min` (default `nothing`).
- id: u_max
  type: Any
  units: n/a
  required: false
  description: Keyword argument `u_max` (default `nothing`).
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
  type: RpoLQMPCController
  units: n/a
  description: Return value of `init_rpo_lqmpc`. Returns `RpoLQMPCController(Ad, Bd,
    horizon, sparse(H), E, F, G, W, model, OSQP.Results()`.
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

# init_rpo_lqmpc

## Purpose
Constructs a ready-to-run RPO LQ-MPC controller from the orbit mean motion, sample time, cost weights and horizon, performing every expensive precomputation once.

## Theory & Math
The horizon cost $J = \tfrac{1}{2}(X - X_{\text{ref}})^\top \bar{Q} (X - X_{\text{ref}}) + \tfrac{1}{2} U^\top \bar{R} U$ with $X = \bar{A} x_0 + \bar{B} U$ reduces to the QP

$$
\min_U \; U^\top H U + 2\,(E x_0 - F X_{\text{ref}})^\top U,\qquad H = \tfrac{1}{2}\left(\bar{B}^\top \bar{Q} \bar{B} + \bar{R}\right)
$$

subject to $u_{\min} \le u_k \le u_{\max}$ for each step $k$.

## Design & Implementation
Chains the continuous CW matrices through zero-order-hold discretisation and horizon stacking, then builds `Qbar` as a block diagonal of `horizon - 1` copies of `Q` followed by `Qf`, and `Rbar` of `horizon` copies of `R`. The Hessian is `0.5 * (Bbar' Qbar Bbar + Rbar)` symmetrised by adding its transpose; `E = Bbar' Qbar Abar` and `F = Bbar' Qbar` give the linear term. Input bounds default to infinite when `u_min` or `u_max` is `nothing` and are expressed as a stacked box constraint through `G` and `W`. OSQP is set up with the upper triangle of `2H`, absolute and relative tolerances of 1e-4, at most 1000 iterations, warm start enabled and polishing off.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n` | Any | n/a | yes | Positional argument `n`. |
| in | `dt` | Any | n/a | yes | Positional argument `dt`. |
| in | `Q` | Any | n/a | yes | Positional argument `Q`. |
| in | `R` | Any | n/a | yes | Positional argument `R`. |
| in | `Qf` | Any | n/a | yes | Positional argument `Qf`. |
| in | `horizon` | Any | n/a | yes | Positional argument `horizon`. |
| in | `u_min` | Any | n/a | no | Keyword argument `u_min` (default `nothing`). |
| in | `u_max` | Any | n/a | no | Keyword argument `u_max` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RpoLQMPCController | n/a | — | Return value of `init_rpo_lqmpc`. Returns `RpoLQMPCController(Ad, Bd, horizon, sparse(H), E, F, G, W, model, OSQP.Results()`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`

**Downstream**

- `callees` → [[core.runtime_types_model|Model]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:105-105`
- `callees` → [[gnc.lqmpc_rpo_block_diag|rpo_block_diag]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:97-97`
- `callees` → [[gnc.lqmpc_rpo_discretize_zoh|rpo_discretize_zoh]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:91-91`
- `callees` → [[gnc.lqmpc_rpo_hcw_continuous_mats|rpo_hcw_continuous_mats]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:90-90`
- `callees` → [[gnc.lqmpc_rpo_prediction_mats|rpo_prediction_mats]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:96-96`
- `callees` → [[gnc.lqmpc_rpolqmpccontroller|RpoLQMPCController]] · `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:121-121`
<!-- vulcan:connections:end -->

## Limitations
Weights are copied to dense `Matrix{Float64}` and the stacking is dense before sparsification, so initialisation memory scales with the square of `nu * horizon`; the solver tolerances and iteration cap are hard-coded, with no keyword to tighten them for a demanding approach.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl` line 89.
