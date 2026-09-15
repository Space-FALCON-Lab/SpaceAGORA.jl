---
id: simulation.dynamics_rhs__all_active_spacecraft_outside_atmosphere
label: _all_active_spacecraft_outside_atmosphere
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _all_active_spacecraft_outside_atmosphere
  lines:
  - 1962
  - 1962
inputs:
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `_all_active_spacecraft_outside_atmosphere`.
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

# _all_active_spacecraft_outside_atmosphere

## Purpose
Lets the implicit-atmosphere RHS skip entirely when every active satellite is above the entry interface under a model that vanishes there.

## Design & Implementation
Loops active satellites and returns false on the first that is not outside the atmosphere per `_spacecraft_outside_atmosphere_for_current_state`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_all_active_spacecraft_outside_atmosphere`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1982-1982`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1965-1965`
<!-- vulcan:connections:end -->

## Limitations
Only models declaring `density_vanishes_above_entry_interface` qualify, so with GRAM this always returns false.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1962.
