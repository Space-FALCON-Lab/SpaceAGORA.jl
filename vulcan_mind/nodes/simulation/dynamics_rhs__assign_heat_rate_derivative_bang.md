---
id: simulation.dynamics_rhs__assign_heat_rate_derivative_bang
label: _assign_heat_rate_derivative!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _assign_heat_rate_derivative!
  lines:
  - 1436
  - 1436
inputs:
- id: du_heat
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `du_heat`.
- id: heat_rates
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `heat_rates`.
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
  description: Return value of `_assign_heat_rate_derivative!`; mutates `du_heat`
    in place. Returns `nothing`.
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

# _assign_heat_rate_derivative!

## Purpose
Copies the per-link heat rates computed by the thermal callback into the heat-load derivative slots of a satellite's state, tolerating a length mismatch between links and slots.

## Design & Implementation
When lengths match it is a single broadcast assignment. Otherwise it zeroes the target and copies the overlapping prefix under `@inbounds`, so extra slots integrate zero and extra rates are dropped. Returns `nothing`. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du_heat` | AbstractVector | n/a | yes | Positional argument `du_heat`. |
| in | `heat_rates` | AbstractVector | n/a | yes | Positional argument `heat_rates`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_assign_heat_rate_derivative!`; mutates `du_heat` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1395-1395`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2147-2147`
- [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1888-1888`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1778-1778`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A mismatch is silently tolerated rather than raised, so a spacecraft whose link count changed after state allocation integrates wrong heat loads without a diagnostic.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1436.
