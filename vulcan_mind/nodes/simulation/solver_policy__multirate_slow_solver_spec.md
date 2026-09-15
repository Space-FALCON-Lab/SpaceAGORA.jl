---
id: simulation.solver_policy__multirate_slow_solver_spec
label: _multirate_slow_solver_spec
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _multirate_slow_solver_spec
  lines:
  - 252
  - 252
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
  description: Return value of `_multirate_slow_solver_spec`. Returns `_multirate_solver_spec_from_sym(cfg.multirate_slow_solver,
    "multirate_slow_solve`.
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

# _multirate_slow_solver_spec

## Purpose
Returns the algorithm spec for the slow (gravity) stage of multirate integration from `cfg.multirate_slow_solver`.

## Design & Implementation
A single-line delegate: `_multirate_solver_spec_from_sym(cfg.multirate_slow_solver, "multirate_slow_solver")`. A zero-argument overload resolves the active `SolverConfig`. The returned NamedTuple has `alg`, `label`, and `auto_switch_capable` fields consumed by `_solve_with_multirate_solver`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | SolverConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_multirate_slow_solver_spec`. Returns `_multirate_solver_spec_from_sym(cfg.multirate_slow_solver, "multirate_slow_solve`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.solver_policy__multirate_fast_solver_spec|_multirate_fast_solver_spec]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:254-254`
- [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:426-426`

**Downstream**

- `callees` → [[simulation.solver_policy__multirate_solver_spec_from_sym|_multirate_solver_spec_from_sym]] · `callers` · call · `src/simulation/engine/solver_policy.jl:252-252`
<!-- vulcan:connections:end -->

## Limitations
Inherits the missing `switch_max` propagation from the underlying mapping. Constructs a fresh algorithm object each call, which is once per run and therefore harmless in practice.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 252.
