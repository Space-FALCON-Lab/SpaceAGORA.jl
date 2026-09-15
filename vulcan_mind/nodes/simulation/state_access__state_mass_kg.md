---
id: simulation.state_access__state_mass_kg
label: _state_mass_kg
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_mass_kg
  lines:
  - 60
  - 60
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  description: Return value of `_state_mass_kg`.
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

# _state_mass_kg

## Purpose
Returns the current mass of spacecraft `sat_idx` in kilograms, bridging the difference between solver modes that propagate mass as a state and the gravity backbone, which does not.

## Design & Implementation
When `_is_gravity_backbone_state(u)` holds, the mass is not part of the propagated state, so the function falls back to the static model and returns `spacecraft.dry_mass + spacecraft.prop_mass` for `args.dynamics_model.spacecraft[sat_idx]`. Otherwise it reads the integrated mass channel `u.sc[sat_idx].mass` and converts it to `Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_state_mass_kg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/state_access.jl:65-65`
- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:61-61`
<!-- vulcan:connections:end -->

## Limitations
In backbone mode the returned value is the full wet mass at model-definition time, so any propellant already expended is ignored and the mass appears frozen for the whole run. This makes backbone results inconsistent with thrust-driven mass depletion. The function also assumes `args.dynamics_model.spacecraft` is indexed identically to the state's `sc` block.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 60.
