---
id: simulation.state_access__state_velocity_ii
label: _state_velocity_ii
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_velocity_ii
  lines:
  - 46
  - 46
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_state_velocity_ii`.
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

# _state_velocity_ii

## Purpose
Reads the inertial velocity vector of spacecraft `sat_idx` out of an integrator state, returning it as an `SVector{3, Float64}` in metres per second regardless of which state layout the active solver uses.

## Design & Implementation
For a gravity-backbone state it pulls the velocity partition via `_gravity_backbone_velocity_state(u)`, unwraps a `.sc` field if present, indexes `sat_idx`, and takes components 1 through 3, because in the second-order layout the velocity block stores velocity in its leading slots. For a flat ComponentVector state it indexes `u.sc[sat_idx]` and takes components 4, 5 and 6, the velocity entries of the packed position-velocity record.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_state_velocity_ii`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`

**Downstream**

- `callees` → [[simulation.state_access__gravity_backbone_velocity_state|_gravity_backbone_velocity_state]] · `callers` · call · `src/simulation/engine/state_access.jl:48-48`
- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
The component offsets are positional constants with no schema check, so any reordering of the flat per-spacecraft state silently returns wrong physics rather than erroring. `sat_idx` is not bounds-checked here, and out-of-range indices surface as an indexing error deeper in ComponentArrays. Values are converted to `Float64` unconditionally, discarding dual numbers if used under forward-mode differentiation.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 46.
