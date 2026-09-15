---
id: simulation.solver_policy__solve_with_explicit_solver
label: _solve_with_explicit_solver
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solve_with_explicit_solver
  lines:
  - 312
  - 312
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
- id: alg
  type: Any
  units: n/a
  required: true
  description: Positional argument `alg`.
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
- id: dtmax_override
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Keyword argument `dtmax_override` (default `nothing`).
- id: solver_cache
  type: Union{Nothing, SolverIntegratorCache}
  units: n/a
  required: false
  description: Keyword argument `solver_cache` (default `nothing`).
- id: needs_full_solution
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `needs_full_solution` (default `true`).
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
  description: Return value of `_solve_with_explicit_solver`. Returns `DiffEqBase.solve!(integ)`
    or `solve(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, save_eve`
    or `solve(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, maxiters`.
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

# _solve_with_explicit_solver

## Purpose
Central adaptive-solve wrapper: resolves `maxiters`, `dtmax`, and save options, then either reuses a cached integrator via `reinit!` or calls `init`/`solve` with the given algorithm and tolerances.

## Design & Implementation
Signature `(prob, cfg, args, alg, reltol_tol, abstol_tol; dtmax_override=nothing, solver_cache=nothing, needs_full_solution=true)`; a second method omits `cfg` and uses `_active_solver_config()`. `dtmax_use` defaults to `dt_max_orbit` and must be positive. `save_everystep` and `save_on` follow the env override when the key exists, else `needs_full_solution`; `save_start`/`save_end` default to `true`. When the cache holds an integrator with matching save flags it sets `integ.p = prob.p`, calls `SciMLBase.reinit!(integ, prob.u0; t0, tf, erase_sol=true, reinit_callbacks=false)`, and returns `DiffEqBase.solve!(integ)`. Otherwise it branches on `maxiters === nothing` and on cache presence, calling `DiffEqBase.init` plus `_cache_integrator!` or plain `solve`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prob` | Any | n/a | yes | Positional argument `prob`. |
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `alg` | Any | n/a | yes | Positional argument `alg`. |
| in | `reltol_tol` | Any | n/a | yes | Positional argument `reltol_tol`. |
| in | `abstol_tol` | Any | n/a | yes | Positional argument `abstol_tol`. |
| in | `dtmax_override` | Union{Nothing, Float64} | n/a | no | Keyword argument `dtmax_override` (default `nothing`). |
| in | `solver_cache` | Union{Nothing, SolverIntegratorCache} | n/a | no | Keyword argument `solver_cache` (default `nothing`). |
| in | `needs_full_solution` | Bool | n/a | no | Keyword argument `needs_full_solution` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_solve_with_explicit_solver`. Returns `DiffEqBase.solve!(integ)` or `solve(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, save_eve` or `solve(prob, alg; reltol=reltol_tol, abstol=abstol_tol, dtmax=dtmax_use, maxiters`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__rodas5p_alg|_rodas5p_alg]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:674-674`
- [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:414-414`
- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:674-674`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/solver_policy.jl:341-341`
- `callees` → [[simulation.from_env__engine_env_haskey_with_env_fallback|_engine_env_haskey_with_env_fallback]] · `callers` · call · `src/simulation/engine/solver_policy.jl:327-327`
- `callees` → [[simulation.solver_policy__cache_integrator_bang|_cache_integrator!]] · `callers` · call · `src/simulation/engine/solver_policy.jl:351-351`
- `callees` → [[simulation.solver_policy__solver_bool_env|_solver_bool_env]] · `callers` · call · `src/simulation/engine/solver_policy.jl:330-330`
- `callees` → [[simulation.solver_policy__solver_cache_options_match|_solver_cache_options_match]] · `callers` · call · `src/simulation/engine/solver_policy.jl:337-337`
- `callees` → [[simulation.solver_policy__solver_save_everystep|_solver_save_everystep]] · `callers` · call · `src/simulation/engine/solver_policy.jl:328-328`
<!-- vulcan:connections:end -->

## Limitations
Cached reuse ignores changes in `alg`, tolerances, `dtmax`, and `maxiters`; the caller must supply a fresh cache when those change. `reinit_callbacks=false` means callback state persists across reinit. The four-way branch duplicates the keyword list, so adding an option requires editing four calls.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 312.
