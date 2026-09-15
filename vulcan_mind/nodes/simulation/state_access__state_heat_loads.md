---
id: simulation.state_access__state_heat_loads
label: _state_heat_loads
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_heat_loads
  lines:
  - 72
  - 72
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
  type: Any
  units: n/a
  description: Return value of `_state_heat_loads`. Returns `zeros(Float64, length(args.dynamics_model.spacecraft[sat_idx].links))`
    or `u.sc[sat_idx].heat_loads`.
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

# _state_heat_loads

## Purpose
Returns the per-link accumulated heat load vector for spacecraft `sat_idx`, used by aerobraking guidance and by thermal reporting.

## Design & Implementation
For a flat ComponentVector state it returns the live `u.sc[sat_idx].heat_loads` array, aliasing solver memory. For a gravity-backbone state, where thermal channels are not propagated, it allocates and returns `zeros(Float64, length(args.dynamics_model.spacecraft[sat_idx].links))` so that the length still matches the spacecraft's link count and downstream indexing stays valid.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_state_heat_loads`. Returns `zeros(Float64, length(args.dynamics_model.spacecraft[sat_idx].links))` or `u.sc[sat_idx].heat_loads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`

**Downstream**

- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:73-73`
<!-- vulcan:connections:end -->

## Limitations
The backbone branch allocates a fresh vector on every call, which is wasteful if invoked inside a derivative or saving loop, and it reports zero heat load rather than signalling that the quantity is unavailable. The flat branch hands back a mutable alias, so a caller that scales or accumulates in place corrupts the integrator state. Length agreement between `links` and the propagated `heat_loads` channel is assumed, never verified.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 72.
