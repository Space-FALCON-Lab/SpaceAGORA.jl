---
id: simulation.execution__build_typed_solver_problem
label: _build_typed_solver_problem
kind: function
source:
  file: src/simulation/engine/execution.jl
  symbol: _build_typed_solver_problem
  lines:
  - 23
  - 23
inputs:
- id: u0
  type: Any
  units: n/a
  required: true
  description: Positional argument `u0`.
- id: tspan
  type: Any
  units: n/a
  required: true
  description: Positional argument `tspan`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: callbacks
  type: Any
  units: n/a
  required: true
  description: Positional argument `callbacks`.
- id: solver_mode
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `solver_mode`.
- id: jac_prototype
  type: Union{Nothing, SparseMatrixCSC{Float64, Int}}
  units: n/a
  required: false
  description: Positional argument `jac_prototype` (default `nothing`).
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
  type: Union{ODEProblem, SecondOrderODEProblem, SplitODEProblem}
  units: n/a
  description: Return value of `_build_typed_solver_problem`. Returns `SecondOrderODEProblem(`
    or `SplitODEProblem(` or `ODEProblem(f, u0, tspan, p; callback=callbacks)` or
    `ODEProblem(spacecraft_dynamics!, u0, tspan, p, callback=callbacks)`.
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

# _build_typed_solver_problem

## Purpose
Selects and constructs the concrete DifferentialEquations.jl problem object that matches the active solver mode, so that `run_simulation` and the checkpoint loop can build a fresh problem for every time segment `(t_cursor, t_next)` with one call.

## Design & Implementation
The five-argument core method switches on `solver_mode::Symbol`. `:gravity_backbone_split` calls `_gravity_backbone_initial_states(u0, p.args)` to obtain `(q0, dq0)` and returns a `SecondOrderODEProblem` over `spacecraft_dynamics_gravity_backbone!`. `:split_imex` returns a `SplitODEProblem` pairing `spacecraft_dynamics_implicit_atmosphere!` with `spacecraft_dynamics_explicit_remainder!`, and `:multirate` pairs `spacecraft_dynamics_slow!` with `spacecraft_dynamics_fast_control!`. Any other mode yields a plain `ODEProblem` over `spacecraft_dynamics!`; when a `jac_prototype::SparseMatrixCSC` is supplied it is wrapped in an `ODEFunction(...; jac_prototype)` first. In all cases `callbacks` is passed as the `callback` keyword. A four-argument convenience method resolves the mode from `_solver_policy_mode(_active_solver_config())`. Both methods are `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u0` | Any | n/a | yes | Positional argument `u0`. |
| in | `tspan` | Any | n/a | yes | Positional argument `tspan`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `callbacks` | Any | n/a | yes | Positional argument `callbacks`. |
| in | `solver_mode` | Symbol | n/a | yes | Positional argument `solver_mode`. |
| in | `jac_prototype` | Union{Nothing, SparseMatrixCSC{Float64, Int}} | n/a | no | Positional argument `jac_prototype` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{ODEProblem, SecondOrderODEProblem, SplitODEProblem} | n/a | — | Return value of `_build_typed_solver_problem`. Returns `SecondOrderODEProblem(` or `SplitODEProblem(` or `ODEProblem(f, u0, tspan, p; callback=callbacks)` or `ODEProblem(spacecraft_dynamics!, u0, tspan, p, callback=callbacks)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:340-340`

**Downstream**

- `callees` → [[simulation.state_access__gravity_backbone_initial_states|_gravity_backbone_initial_states]] · `callers` · call · `src/simulation/engine/execution.jl:27-27`
<!-- vulcan:connections:end -->

## Limitations
Unknown mode symbols are not rejected; they fall through to the default `ODEProblem`, so a typo in `solver_mode` silently changes the integration scheme. The `jac_prototype` is ignored for the three split modes even if provided. The return type varies with mode (three different problem types), which forces type instability at the call site and is why the checkpoint loop rebuilds the problem per segment rather than reusing it.

## Provenance
Mapped from `src/simulation/engine/execution.jl` line 23.
