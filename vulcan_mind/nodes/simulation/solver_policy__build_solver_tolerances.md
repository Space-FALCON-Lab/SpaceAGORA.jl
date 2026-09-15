---
id: simulation.solver_policy__build_solver_tolerances
label: _build_solver_tolerances
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _build_solver_tolerances
  lines:
  - 15
  - 15
inputs:
- id: u_state
  type: ComponentVector
  units: n/a
  required: true
  description: Positional argument `u_state`.
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
  type: Any
  units: n/a
  description: Return value of `_build_solver_tolerances`. Returns `tol.reltol_orbit,
    tol.abstol_orbit` or `reltol_state, abstol_state`.
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

# _build_solver_tolerances

## Purpose
Produces either scalar `(reltol, abstol)` or full `ComponentVector` tolerance arrays matching the state layout `u_state`, assigning distinct tolerances to mass, heat-load, and attitude components.

## Design & Implementation
Returns `(tol.reltol_orbit, tol.abstol_orbit)` when `_requires_componentwise_tolerances(args)` is false. Otherwise it resolves six component tolerances via `_resolve_component_tolerance`, copies `u_state` twice, fills both with the orbit values, and in one `@inbounds` loop over `reltol_state.sc` sets `.mass`, `.heat_loads .=`, and, when `orientation_sim`, `.ω .=` and `.q .=` (quaternion tolerances come straight from `tol.reltol_quaternion`/`abstol_quaternion` without the zero-fallback). Returns `(reltol_state, abstol_state)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u_state` | ComponentVector | n/a | yes | Positional argument `u_state`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_build_solver_tolerances`. Returns `tol.reltol_orbit, tol.abstol_orbit` or `reltol_state, abstol_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:302-302`

**Downstream**

- `callees` → [[simulation.solver_policy__requires_componentwise_tolerances|_requires_componentwise_tolerances]] · `callers` · call · `src/simulation/engine/solver_policy.jl:17-17`
- `callees` → [[simulation.solver_policy__resolve_component_tolerance|_resolve_component_tolerance]] · `callers` · call · `src/simulation/engine/solver_policy.jl:21-21`
<!-- vulcan:connections:end -->

## Limitations
Two full copies of the state vector are allocated per call. Quaternion tolerances bypass `_resolve_component_tolerance`, so a `0.0` quaternion tolerance is passed literally to the solver. The function assumes each `sc[i]` block has `mass`, `heat_loads`, `ω`, and `q` fields; a layout without `heat_loads` throws on the broadcast.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 15.
