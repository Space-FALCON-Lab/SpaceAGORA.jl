---
id: simulation.state_access__state_has_mass
label: _state_has_mass
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_has_mass
  lines:
  - 68
  - 68
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
  description: Return value of `_state_has_mass`.
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

# _state_has_mass

## Purpose
Capability predicate answering whether a mass value can be obtained for spacecraft `sat_idx` from the given state and model arguments.

## Design & Implementation
Unconditionally returns `true`, ignoring all three of its arguments `u`, `args` and `sat_idx`. This is correct by construction because its partner `_state_mass_kg` always produces a number: in the flat layout it reads the integrated `mass` channel, and in the gravity-backbone layout it falls back to the static sum of `dry_mass` and `prop_mass` from the dynamics model. The predicate exists so that saving and reporting code can use a uniform has-value/get-value pair across state fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_state_has_mass`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the answer is constant, it conveys no information and cannot warn callers that the backbone path is returning a static model mass rather than a propagated one. Consumers that treat a true result as evidence of an integrated mass state will silently record a constant mass timeline.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 68.
