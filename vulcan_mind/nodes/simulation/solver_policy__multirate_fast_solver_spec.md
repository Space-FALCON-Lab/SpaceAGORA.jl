---
id: simulation.solver_policy__multirate_fast_solver_spec
label: _multirate_fast_solver_spec
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _multirate_fast_solver_spec
  lines:
  - 253
  - 253
inputs:
- id: cfg
  type: SolverConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_multirate_fast_solver_spec`. Returns `_multirate_solver_spec_from_sym(cfg.multirate_fast_solver,
    "multirate_fast_solve`.
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

# _multirate_fast_solver_spec

## Purpose
Returns the algorithm spec for the fast (perturbation) stage of multirate integration from `cfg.multirate_fast_solver`.

## Design & Implementation
Delegates to `_multirate_solver_spec_from_sym(cfg.multirate_fast_solver, "multirate_fast_solver")`; the zero-argument overload reads `_active_solver_config()`. The `alg` is used for both fast half-steps of each Strang macro step with `dtmax_override=min(fast_dt_s, half_dt)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_multirate_fast_solver_spec`. Returns `_multirate_solver_spec_from_sym(cfg.multirate_fast_solver, "multirate_fast_solve`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:427-427`

**Downstream**

- `callees` → [[simulation.solver_policy__multirate_slow_solver_spec|_multirate_slow_solver_spec]] · `callers` · call · `src/simulation/engine/solver_policy.jl:254-254`
- `callees` → [[simulation.solver_policy__multirate_solver_spec_from_sym|_multirate_solver_spec_from_sym]] · `callers` · call · `src/simulation/engine/solver_policy.jl:253-253`
<!-- vulcan:connections:end -->

## Limitations
Same restrictions as the slow spec: the `:auto_stiff` choice uses the library default `switch_max`. Because the fast stage is invoked twice per macro step, an implicit choice here doubles the Jacobian cost per macro step.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 253.
