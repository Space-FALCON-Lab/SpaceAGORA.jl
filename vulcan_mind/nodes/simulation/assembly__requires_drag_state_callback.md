---
id: simulation.assembly__requires_drag_state_callback
label: _requires_drag_state_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_drag_state_callback
  lines:
  - 73
  - 73
inputs:
- id: effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `effectors`.
- id: args
  type: SimulationConfiguration
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
  type: Bool
  units: n/a
  description: Return value of `_requires_drag_state_callback`.
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

# _requires_drag_state_callback

## Purpose
Determines whether the run needs the drag-state switching callback, which swaps the integrator's tolerance and step-size limits between atmospheric and orbital flight phases.

## Design & Implementation
Returns `false` immediately when `_requires_density_callback(effectors, args)` fails. Otherwise it binds `tol = args.integration_tolerances` and returns `true` when any of the three atmosphere-versus-orbit pairs differ: `tol.dt_max_atmosphere != tol.dt_max_orbit`, `tol.reltol_atmosphere != tol.reltol_orbit`, or `tol.abstol_atmosphere != tol.abstol_orbit`. The reasoning is that if all three pairs are identical there is nothing to switch, so the callback and its staged-density dependency can both be dropped.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_drag_state_callback`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.assembly__requires_staged_density_callback|_requires_staged_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:40-40`
- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:176-176`

**Downstream**

- `callees` → [[simulation.assembly__requires_density_callback|_requires_density_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:74-74`
<!-- vulcan:connections:end -->

## Limitations
The comparisons use exact floating-point inequality, so tolerances that differ only in the last bit — for example from a round-trip through a text scenario file — count as different and install the callback for no practical benefit, while tolerances intended to differ but accidentally written identically silently disable phase switching. Only these three pairs are checked: a scenario that differentiates a component tolerance such as `reltol_mass` between phases, but leaves the three baseline pairs equal, gets no switching. The predicate is structural and does not check that the atmospheric tolerances are actually tighter than the orbital ones.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 73.
