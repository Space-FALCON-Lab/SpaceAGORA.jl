---
id: simulation.assembly__resolved_component_tolerance
label: _resolved_component_tolerance
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _resolved_component_tolerance
  lines:
  - 95
  - 95
inputs:
- id: component_tol
  type: Float64
  units: n/a
  required: true
  description: Positional argument `component_tol`.
- id: baseline_tol
  type: Float64
  units: n/a
  required: true
  description: Positional argument `baseline_tol`.
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
  description: Return value of `_resolved_component_tolerance`.
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

# _resolved_component_tolerance

## Purpose
Applies the unset-means-inherit convention for per-component integration tolerances, substituting the phase baseline whenever a component tolerance was left at its zero default.

## Design & Implementation
Written `@inline _resolved_component_tolerance(component_tol, baseline_tol) = component_tol == 0.0 ? baseline_tol : component_tol`, with both arguments and the result typed `Float64`. `_callback_tolerances_for_phase` calls it six times, for the mass, heat-load, and angular-rate relative and absolute tolerances, against whichever of the atmospheric or orbital baselines the current phase selects.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `component_tol` | Float64 | n/a | yes | Positional argument `component_tol`. |
| in | `baseline_tol` | Float64 | n/a | yes | Positional argument `baseline_tol`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_resolved_component_tolerance`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- [[simulation.assembly__callback_tolerances_for_phase|_callback_tolerances_for_phase]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:107-107`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Zero is overloaded as the sentinel for unset, so a user who deliberately wants a zero absolute tolerance on one component — a legitimate request meaning pure relative control — cannot express it: the value is silently replaced by the baseline. The comparison `component_tol == 0.0` also matches negative zero, and a negative tolerance is passed through unchanged rather than rejected. No validation ensures the resolved value is finite.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 95.
