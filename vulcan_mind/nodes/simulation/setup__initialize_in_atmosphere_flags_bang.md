---
id: simulation.setup__initialize_in_atmosphere_flags_bang
label: _initialize_in_atmosphere_flags!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_in_atmosphere_flags!
  lines:
  - 1343
  - 1343
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: initial_conditions
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_conditions`.
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
  type: Nothing
  units: n/a
  description: Return value of `_initialize_in_atmosphere_flags!`; mutates `p` in
    place.
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

# _initialize_in_atmosphere_flags!

## Purpose
Seeds each satellite's `in_atmosphere` flag from its initial altitude so the first RHS evaluation has a correct value before any callback has fired.

## Design & Implementation
Resizes the flag vector to the satellite count, computes altitude as the position norm minus `Rp_e`, and sets the flag to `alt <= EI` in metres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `initial_conditions` | Any | n/a | yes | Positional argument `initial_conditions`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_in_atmosphere_flags!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:185-185`

**Downstream**

- `callees` → [[simx.engine_state_access_state_position_ii|_state_position_ii]] · `callers` · call · `src/simulation/engine/setup.jl:1350-1350`
<!-- vulcan:connections:end -->

## Limitations
Uses spherical altitude for this initial test while the callbacks later use geodetic altitude, so a satellite starting within the flattening band of the interface can flip on the first callback.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1343.
