---
id: simulation.solver_policy__symplectic_fixed_dt_s
label: _symplectic_fixed_dt_s
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _symplectic_fixed_dt_s
  lines:
  - 86
  - 86
inputs:
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
  type: Float64
  units: n/a
  description: Return value of `_symplectic_fixed_dt_s`.
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

# _symplectic_fixed_dt_s

## Purpose
Returns the fixed step (seconds) for the symplectic `KahanLi8` mode, defaulting to `dt_max_orbit` when `SolverConfig.symplectic_dt_s` is `nothing`.

## Design & Implementation
`dt = isnothing(cfg.symplectic_dt_s) ? args.integration_tolerances.dt_max_orbit : cfg.symplectic_dt_s`; throws `ArgumentError("SolverConfig.symplectic_dt_s must be > 0.0, got ...")` unless `dt > 0.0`. A one-argument overload resolves `cfg` via `_active_solver_config()`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_symplectic_fixed_dt_s`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:642-642`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Falling back to `dt_max_orbit` (a cap intended for adaptive solvers) can yield a very coarse symplectic step with no warning. NaN passes the `> 0.0` guard as false and is rejected, but `Inf` passes and produces a single giant step. No check that `dt` divides the span.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 86.
