---
id: simulation.save_fields__save_lift
label: _save_lift
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_lift
  lines:
  - 43
  - 43
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `_save_lift`. Returns `lifts`.
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

# _save_lift

## Purpose
Save-time getter for the aerodynamic lift force on each spacecraft, taken from the lift cache filled by the aerodynamic effector during right-hand-side evaluation.

## Design & Implementation
Marked `@inline` and identical in shape to the drag getter, but sourcing `integrator.p.save_cache.lift_cache`. It builds a `Vector{SVector{3, Float64}}` of length `num_sats`, copying `lift_cache[i]` where the index is in range and writing `SVector{3, Float64}(0.0, 0.0, 0.0)` otherwise. Lift is cached separately from drag so that the saved decomposition matches the effector's own along-flow and normal split rather than being reconstructed from a total force.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_save_lift`. Returns `lifts`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:181-181`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The same staleness applies as for drag: the cached vector is from the last right-hand-side call, not necessarily the accepted step. A spacecraft configured with an aerodynamic model that does not populate the lift cache saves as exactly zero with no indication that the quantity was never computed. The frame of the cached vector is set by the effector and is not restated or converted here.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 43.
