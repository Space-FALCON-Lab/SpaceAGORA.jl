---
id: simulation.assembly__requires_density_for_rhs
label: _requires_density_for_rhs
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _requires_density_for_rhs
  lines:
  - 19
  - 19
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
  description: Return value of `_requires_density_for_rhs`.
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

# _requires_density_for_rhs

## Purpose
The canonical test for whether atmospheric density is needed anywhere in the dynamics right-hand side, combining the presence of an aerodynamic effector with the presence of a non-trivial atmosphere model.

## Design & Implementation
Returns `_uses_atmospheric_dynamic_effector(effectors) || !(args.environment_model.density_model isa NoAtmosphereModel)`. The disjunction is deliberately generous: an aerodynamic effector implies density is consumed, while any density model other than the `NoAtmosphereModel` sentinel implies density is produced and may be consumed by a thermal or entry-detection path even without an aero force. Every other requirement predicate in this file, including `_requires_thermal_callback`, `_requires_drag_state_callback`, and `_requires_entry_end_callback`, short-circuits on this result first.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_requires_density_for_rhs`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- [[simulation.assembly__requires_density_callback|_requires_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:26-26`
- [[simulation.assembly__requires_staged_density_callback|_requires_staged_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:34-34`

**Downstream**

- `callees` → [[simulation.assembly__uses_atmospheric_dynamic_effector|_uses_atmospheric_dynamic_effector]] · `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:20-20`
<!-- vulcan:connections:end -->

## Limitations
Because it is an `or`, a scenario that configures a real GRAM atmosphere but no aerodynamic effector still reports `true` and pays for density machinery it never reads. Conversely, the aero-effector test inherits the hard-coded type list of `_uses_atmospheric_dynamic_effector`, so an unrecognised aero model combined with `NoAtmosphereModel` yields `false` and the run proceeds with density unavailable, failing later inside the effector. The predicate is purely structural and does not consider altitude, so a purely exo-atmospheric orbit still triggers the full atmospheric path.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 19.
