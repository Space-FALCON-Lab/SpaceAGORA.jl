---
id: gnc.lqmpc_rpolqmpccontroller
label: RpoLQMPCController
kind: struct
source:
  file: src/gnc/control/rpo_mpc/lqmpc.jl
  symbol: RpoLQMPCController
  lines:
  - 72
  - 72
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
  type: SparseMatrixCSC{Float64, Int}
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
- id: G
  type: SparseMatrixCSC{Float64, Int}
  units: n/a
  required: true
  description: Field `G`.
- id: W
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `W`.
- id: qp_model
  type: OSQP.Model
  units: n/a
  required: true
  description: Field `qp_model`.
- id: qp_results
  type: OSQP.Results
  units: n/a
  required: true
  description: Field `qp_results`.
- id: U_prev
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `U_prev`.
- id: u_min
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `u_min`.
- id: u_max
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `u_max`.
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
  description: Constructed `RpoLQMPCController`.
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

# RpoLQMPCController

## Purpose
Holds every precomputed quantity the RPO tracking MPC needs per tick — discrete dynamics, cost matrices, input bounds, the live OSQP model and the warm-start vector — so the online step is one QP update and solve.

## Design & Implementation
A `mutable struct` of thirteen fields. `Ad`, `Bd` and `horizon` describe the prediction model; `H` is the sparse Hessian, `E` and `F` the dense linear-term maps from initial state and stacked reference, and `G` and `W` the sparse box-constraint matrix and bound vector. `qp_model` and `qp_results` are the OSQP handles reused across solves, `U_prev` is the shifted previous solution used as the warm start, and `u_min` and `u_max` retain the per-axis bounds. Mutability is needed only for `U_prev`, which the control step overwrites in place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `Ad` | Matrix{Float64} | n/a | yes | Field `Ad`. |
| in | `Bd` | Matrix{Float64} | n/a | yes | Field `Bd`. |
| in | `horizon` | Int | n/a | yes | Field `horizon`. |
| in | `H` | SparseMatrixCSC{Float64, Int} | n/a | yes | Field `H`. |
| in | `E` | Matrix{Float64} | n/a | yes | Field `E`. |
| in | `F` | Matrix{Float64} | n/a | yes | Field `F`. |
| in | `G` | SparseMatrixCSC{Float64, Int} | n/a | yes | Field `G`. |
| in | `W` | Vector{Float64} | n/a | yes | Field `W`. |
| in | `qp_model` | OSQP.Model | n/a | yes | Field `qp_model`. |
| in | `qp_results` | OSQP.Results | n/a | yes | Field `qp_results`. |
| in | `U_prev` | Vector{Float64} | n/a | yes | Field `U_prev`. |
| in | `u_min` | Vector{Float64} | n/a | yes | Field `u_min`. |
| in | `u_max` | Vector{Float64} | n/a | yes | Field `u_max`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RpoLQMPCController | n/a | — | Constructed `RpoLQMPCController`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.lqmpc_init_rpo_lqmpc|init_rpo_lqmpc]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:121-121`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct owns a native OSQP model, so it cannot be serialized or deep-copied safely, and two tasks sharing one controller would race on `qp_model` and `U_prev`; the bounds are stored but the constraint set is fixed at setup, so changing `u_min` after construction has no effect on the solver.

## Provenance
Mapped from `src/gnc/control/rpo_mpc/lqmpc.jl` line 72.
