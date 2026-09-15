---
id: simulation.solver_policy__solve_with_multirate_solver
label: _solve_with_multirate_solver
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solve_with_multirate_solver
  lines:
  - 406
  - 406
inputs:
- id: prob
  type: Any
  units: n/a
  required: true
  description: Positional argument `prob`.
- id: cfg
  type: SolverConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: reltol_tol
  type: Any
  units: n/a
  required: true
  description: Positional argument `reltol_tol`.
- id: abstol_tol
  type: Any
  units: n/a
  required: true
  description: Positional argument `abstol_tol`.
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
  description: Return value of `_solve_with_multirate_solver`. Returns `sol, (` or
    `sol_fast_pre, (` or `sol_slow, (` or `sol_fast_post, (` (and more).
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

# _solve_with_multirate_solver

## Purpose
Integrates a split problem with second-order Strang splitting: fast half-step on `f2`, slow full step on `f1`, fast half-step on `f2`, repeated over macro steps of `slow_dt_s`, returning the last segment solution and metadata.

## Theory & Math
One macro step of length $h$ applies $u_{n+1} = \Phi^{f_2}_{h/2} \circ \Phi^{f_1}_{h} \circ \Phi^{f_2}_{h/2}(u_n)$, where $\Phi^{f}_{\tau}$ denotes the exact-in-time flow of $\dot u = f(u,t)$ over $\tau$ seconds approximated by the configured sub-solver; the composition is second-order accurate in $h$.

## Design & Implementation
Requires `prob.f` to have `f1` and `f2`, else throws `ArgumentError`. For an empty span it falls back to a single `Tsit5` solve. Otherwise it resolves slow and fast specs, `fast_substeps`, `slow_dt_s`, and `fast_dt_s = slow_dt_s / fast_substeps`. The loop advances `t_cursor` by `min(slow_dt_s, remaining)`, building subproblems with `_split_subproblem` and solving with `_solve_with_explicit_solver` using `dtmax_override=min(fast_dt_s, half_dt)` for fast stages and `segment_dt` for the slow stage. After each stage it counts `_auto_stiff_switched` events for auto-capable specs, returns early on an unsuccessful retcode, and `deepcopy`s `sol.u[end]` into `u_cursor`. Returns `(final_sol, (slow_solver, fast_solver, macro_steps, fast_substeps, slow_dt_s, fast_dt_s, auto_switch_events))`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prob` | Any | n/a | yes | Positional argument `prob`. |
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `reltol_tol` | Any | n/a | yes | Positional argument `reltol_tol`. |
| in | `abstol_tol` | Any | n/a | yes | Positional argument `abstol_tol`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_solve_with_multirate_solver`. Returns `sol, (` or `sol_fast_pre, (` or `sol_slow, (` or `sol_fast_post, (` (and more). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__rodas5p_alg|_rodas5p_alg]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:721-721`
- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:721-721`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/solver_policy.jl:411-411`
- `callees` → [[simulation.solver_policy__auto_stiff_switched|_auto_stiff_switched]] · `callers` · call · `src/simulation/engine/solver_policy.jl:458-458`
- `callees` → [[simulation.solver_policy__multirate_fast_solver_spec|_multirate_fast_solver_spec]] · `callers` · call · `src/simulation/engine/solver_policy.jl:427-427`
- `callees` → [[simulation.solver_policy__multirate_slow_dt_s|_multirate_slow_dt_s]] · `callers` · call · `src/simulation/engine/solver_policy.jl:429-429`
- `callees` → [[simulation.solver_policy__multirate_slow_solver_spec|_multirate_slow_solver_spec]] · `callers` · call · `src/simulation/engine/solver_policy.jl:426-426`
- `callees` → [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callers` · call · `src/simulation/engine/solver_policy.jl:414-414`
- `callees` → [[simulation.solver_policy__split_subproblem|_split_subproblem]] · `callers` · call · `src/simulation/engine/solver_policy.jl:448-448`
<!-- vulcan:connections:end -->

## Limitations
Only the final sub-segment's solution is returned, so trajectory output covers at most the last macro step; callers needing a full history cannot use this mode. Three `deepcopy` calls per macro step plus three problem constructions allocate heavily. `fast_substeps` is not validated for positivity, so `0` yields `fast_dt_s = Inf`. No `solver_cache` is threaded through, so every stage re-initialises its integrator.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 406.
