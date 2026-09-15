---
id: simulation.execution__build_block_diagonal_jac_prototype
label: _build_block_diagonal_jac_prototype
kind: function
source:
  file: src/simulation/engine/execution.jl
  symbol: _build_block_diagonal_jac_prototype
  lines:
  - 1
  - 1
inputs:
- id: u0
  type: ComponentVector
  units: n/a
  required: true
  description: Positional argument `u0`.
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
  type: SparseMatrixCSC{Float64,
  units: n/a
  description: Return value of `_build_block_diagonal_jac_prototype`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _build_block_diagonal_jac_prototype

## Purpose
Constructs the sparsity pattern of the Jacobian for a multi-satellite state vector so that implicit ODE solvers can exploit the fact that satellites in an uncoupled constellation do not influence each other's derivatives. The result is a `SparseMatrixCSC{Float64,Int}` with one dense block per spacecraft on the diagonal and zeros elsewhere.

## Design & Implementation
Takes the initial `ComponentVector` `u0` and reads `n_sats = length(u0.sc)` and `n_total = length(u0)`. For each satellite `i` it builds `probe = zero(u0)`, sets `probe.sc[i] .= 1.0`, flattens it to a `Vector{Float64}` and uses `findall(!iszero, ...)` to recover the flat index set of that satellite's block. Every ordered pair `(k, j)` from that set is pushed to `rows`/`cols`, and `sparse(rows, cols, ones(...), n_total, n_total)` assembles the matrix. Probing rather than arithmetic on block sizes lets satellites carry heterogeneous state lengths (different `n_bodies` or heat-load counts). `run_simulation` calls it only when `solver_mode` is not one of `:gravity_backbone_split`, `:split_imex`, `:multirate` and `length(u_start.sc) > 1`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u0` | ComponentVector | n/a | yes | Positional argument `u0`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SparseMatrixCSC{Float64, | n/a | — | Return value of `_build_block_diagonal_jac_prototype`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:308-308`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/execution.jl:16-16`
<!-- vulcan:connections:end -->

## Limitations
Each block costs `n_i^2` index pushes, so memory grows quadratically with per-satellite state size; a large multibody satellite produces a dense block regardless of its true internal sparsity. Any state components outside `u0.sc` (shared or global states) are left out of the pattern entirely, which would make a Newton iteration treat their couplings as zero. The pattern also assumes satellites are truly uncoupled; dynamic effectors that couple satellites violate it silently.

## Provenance
Mapped from `src/simulation/engine/execution.jl` line 1.
