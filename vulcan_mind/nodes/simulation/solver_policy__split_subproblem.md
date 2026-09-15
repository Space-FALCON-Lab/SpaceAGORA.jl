---
id: simulation.solver_policy__split_subproblem
label: _split_subproblem
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _split_subproblem
  lines:
  - 402
  - 402
inputs:
- id: prob
  type: Any
  units: n/a
  required: true
  description: Positional argument `prob`.
- id: f
  type: Any
  units: n/a
  required: true
  description: Positional argument `f`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: tspan
  type: Any
  units: n/a
  required: true
  description: Positional argument `tspan`.
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
  type: ODEProblem
  units: n/a
  description: Return value of `_split_subproblem`. Returns `ODEProblem(f, u, tspan,
    prob.p; prob.kwargs...)`.
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

# _split_subproblem

## Purpose
Builds a plain `ODEProblem` for one half of a split RHS (`f1` or `f2`) over a sub-interval, reusing the parent problem's parameters and keyword arguments.

## Design & Implementation
Returns `ODEProblem(f, u, tspan, prob.p; prob.kwargs...)`. Called three times per macro step in `_solve_with_multirate_solver` with `prob.f.f2` for the fast half-steps and `prob.f.f1` for the slow full step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `prob` | Any | n/a | yes | Positional argument `prob`. |
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `tspan` | Any | n/a | yes | Positional argument `tspan`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ODEProblem | n/a | — | Return value of `_split_subproblem`. Returns `ODEProblem(f, u, tspan, prob.p; prob.kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:448-448`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Forwarding `prob.kwargs` re-attaches any callbacks to every subproblem, so discrete callbacks fire on each sub-interval boundary. Allocates a new problem object per stage. The `u` passed is the caller's `u_cursor`, which the caller must `deepcopy` to avoid aliasing.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 402.
