---
id: simulation.state_access__state_has_heat_loads
label: _state_has_heat_loads
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_has_heat_loads
  lines:
  - 79
  - 79
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
  type: Bool
  units: n/a
  description: Return value of `_state_has_heat_loads`.
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

# _state_has_heat_loads

## Purpose
Predicate reporting whether the integrator state actually carries propagated per-link heat loads for spacecraft `sat_idx`, so that saving code can skip thermal columns when the backbone solver is running.

## Design & Implementation
Implemented as `!_is_gravity_backbone_state(u)`: the flat ComponentVector layout includes a `heat_loads` channel per spacecraft, while the two-partition second-order backbone layout carries only position and velocity. The `args` and `sat_idx` arguments are accepted for signature symmetry with `_state_heat_loads` but are not consulted.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_state_has_heat_loads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`

**Downstream**

- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:80-80`
<!-- vulcan:connections:end -->

## Limitations
The check is at the whole-state level, so it cannot express the case where one spacecraft in a heterogeneous fleet lacks thermal channels while another has them. It also does not verify that the `heat_loads` field is actually present on a flat state, so a reduced flat state without thermal channels would be reported as having heat loads and fail later at the property access.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 79.
