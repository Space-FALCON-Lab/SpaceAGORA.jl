---
id: simulation.assembly__callback_tolerances_for_phase
label: _callback_tolerances_for_phase
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _callback_tolerances_for_phase
  lines:
  - 99
  - 99
inputs:
- id: template_reltol
  type: Any
  units: n/a
  required: true
  description: Positional argument `template_reltol`.
- id: template_abstol
  type: Any
  units: n/a
  required: true
  description: Positional argument `template_abstol`.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: in_atmosphere
  type: Bool
  units: n/a
  required: true
  description: Positional argument `in_atmosphere`.
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
  description: Return value of `_callback_tolerances_for_phase`. Returns `baseline_reltol,
    baseline_abstol` or `reltol_new, abstol_new`.
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

# _callback_tolerances_for_phase

## Purpose
Builds the integrator relative- and absolute-tolerance objects for the current flight phase, selecting atmospheric or orbital baselines and then overriding per-component tolerances for mass, heat loads, quaternions, and angular rates.

## Design & Implementation
Picks `baseline_reltol` and `baseline_abstol` from `tol.reltol_atmosphere`/`tol.abstol_atmosphere` when `in_atmosphere` is true and the `*_orbit` fields otherwise. If both templates are plain `Number`s there is no structure to fill, so the two scalars are returned immediately. Otherwise the six component values are resolved through `_resolved_component_tolerance`, the templates are `copy`d, and the copies are broadcast-filled with the baselines using `.=` before per-field overrides are written. The `@inbounds` loop over `eachindex(reltol_new.sc)` sets `sc[i].mass` scalar-wise, broadcasts into `sc[i].heat_loads`, and uses `hasproperty` guards to write `sc[i].q` from `tol.reltol_quaternion`/`tol.abstol_quaternion` and `sc[i].ω` from the resolved angular-rate values only when those fields exist. The `hasproperty` guards are what let one code path serve both attitude-enabled and translation-only state layouts.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `template_reltol` | Any | n/a | yes | Positional argument `template_reltol`. |
| in | `template_abstol` | Any | n/a | yes | Positional argument `template_abstol`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `in_atmosphere` | Bool | n/a | yes | Positional argument `in_atmosphere`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_callback_tolerances_for_phase`. Returns `baseline_reltol, baseline_abstol` or `reltol_new, abstol_new`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- [[simulation.event_callbacks_affect_upcrossing_bang|affect_upcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:161-161`

**Downstream**

- `callees` → [[simulation.assembly__resolved_component_tolerance|_resolved_component_tolerance]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:107-107`
<!-- vulcan:connections:end -->

## Limitations
`copy` on the template is shallow, so if the tolerance state object nests arrays that `copy` does not duplicate, the broadcast writes mutate the caller's template in place — a real hazard given the function is called repeatedly with the same template. Quaternion tolerances bypass `_resolved_component_tolerance` entirely and are written raw, so a zero there is taken literally while a zero mass tolerance inherits the baseline: the two conventions are inconsistent. The `hasproperty` checks are inside the per-spacecraft loop and repeat for every vehicle even though the layout is uniform. Nothing validates that the resolved tolerances are positive, so a misconfigured zero relative tolerance reaches the solver.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 99.
