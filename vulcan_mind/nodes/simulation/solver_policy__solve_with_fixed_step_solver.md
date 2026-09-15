---
id: simulation.solver_policy__solve_with_fixed_step_solver
label: _solve_with_fixed_step_solver
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _solve_with_fixed_step_solver
  lines:
  - 381
  - 381
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
- id: alg
  type: Any
  units: n/a
  required: true
  description: Positional argument `alg`.
- id: dt_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt_s`.
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
  description: Return value of `_solve_with_fixed_step_solver`. Returns `solve(prob,
    alg; dt=dt_s, save_everystep=save_everystep, save_on=save_on, save_s` or `solve(prob,
    alg; dt=dt_s, maxiters=maxiters, save_everystep=save_everystep, save`.
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

# _solve_with_fixed_step_solver

## Purpose
Runs a fixed-step solve (used for `KahanLi8` symplectic and gravity-backbone cores) with the resolved save options and optional `maxiters`.

## Design & Implementation
Resolves `maxiters` from `cfg`, `save_everystep`/`save_on` from env override or `needs_full_solution`, and `save_start`/`save_end` from `_solver_bool_env`. Calls `solve(prob, alg; dt=dt_s, ...)`, adding `maxiters` only when it is not `nothing`. A two-argument convenience overload `(prob, alg, dt_s)` uses the active config and default `needs_full_solution=true`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prob` | Any | n/a | yes | Positional argument `prob`. |
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `alg` | Any | n/a | yes | Positional argument `alg`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `needs_full_solution` | Bool | n/a | no | Keyword argument `needs_full_solution` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_solve_with_fixed_step_solver`. Returns `solve(prob, alg; dt=dt_s, save_everystep=save_everystep, save_on=save_on, save_s` or `solve(prob, alg; dt=dt_s, maxiters=maxiters, save_everystep=save_everystep, save`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:585-585`
- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:642-642`

**Downstream**

- `callees` → [[simulation.from_env__engine_env_haskey_with_env_fallback|_engine_env_haskey_with_env_fallback]] · `callers` · call · `src/simulation/engine/solver_policy.jl:384-384`
- `callees` → [[simulation.solver_policy__solver_bool_env|_solver_bool_env]] · `callers` · call · `src/simulation/engine/solver_policy.jl:387-387`
- `callees` → [[simulation.solver_policy__solver_save_everystep|_solver_save_everystep]] · `callers` · call · `src/simulation/engine/solver_policy.jl:385-385`
<!-- vulcan:connections:end -->

## Limitations
No tolerance arguments are accepted, so an adaptive `alg` passed here would run with library defaults. `dt_s` is not validated for positivity in this function; callers rely on `_symplectic_fixed_dt_s` or `_gravity_backbone_fixed_dt_s`. The convenience overload cannot disable full-solution saving.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 381.
