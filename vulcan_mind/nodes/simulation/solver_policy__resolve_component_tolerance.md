---
id: simulation.solver_policy__resolve_component_tolerance
label: _resolve_component_tolerance
kind: function
source:
  file: src/simulation/engine/solver_policy.jl
  symbol: _resolve_component_tolerance
  lines:
  - 1
  - 1
inputs:
- id: component_tol
  type: Float64
  units: n/a
  required: true
  description: Positional argument `component_tol`.
- id: fallback_tol
  type: Float64
  units: n/a
  required: true
  description: Positional argument `fallback_tol`.
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
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
  description: Return value of `_resolve_component_tolerance`.
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

# _resolve_component_tolerance

## Purpose
Resolves one per-component integration tolerance by treating `0.0` as "inherit the orbit tolerance" and rejecting negative values.

## Design & Implementation
Takes `component_tol::Float64`, `fallback_tol::Float64`, and a `name::String` used only for the error message. Throws `ArgumentError("$name must be >= 0.0, got ...")` when `component_tol < 0.0`; returns `fallback_tol` when the component value is exactly `0.0`; otherwise returns `component_tol`. Marked `@inline` and called six times from `_build_solver_tolerances` for mass, heat load, and angular-rate reltol/abstol.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `component_tol` | Float64 | n/a | yes | Positional argument `component_tol`. |
| in | `fallback_tol` | Float64 | n/a | yes | Positional argument `fallback_tol`. |
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_resolve_component_tolerance`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/solver_policy.jl`
- [[simulation.solver_policy__build_solver_tolerances|_build_solver_tolerances]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:21-21`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Using `0.0` as a sentinel means a genuinely zero absolute tolerance cannot be requested. The fallback is not validated, so a negative `reltol_orbit` passes through. NaN inputs fail the `< 0.0` test and the `== 0.0` test, so NaN is returned as a tolerance.

## Provenance
Mapped from `src/simulation/engine/solver_policy.jl` line 1.
