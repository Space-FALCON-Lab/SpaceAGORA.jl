---
id: simulation.save_fields__save_wind
label: _save_wind
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_wind
  lines:
  - 63
  - 63
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
  description: Return value of `_save_wind`. Returns `winds`.
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

# _save_wind

## Purpose
Save-time getter for the atmospheric wind vector sampled at each spacecraft, recording what the atmosphere model actually supplied to the aerodynamic computation.

## Design & Implementation
Marked `@inline`. Unlike the force getters it reads from `integrator.p.shared_buffers.winds` rather than `save_cache`, because wind is written by the density callbacks into the shared atmosphere buffers and consumed by both aerodynamics and the thermal path. It returns a `Vector{SVector{3, Float64}}` of length `num_sats`, guarding each read with `i <= length(shared_winds)` and substituting the zero vector for out-of-range indices.

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
| out | `result` | Any | n/a | — | Return value of `_save_wind`. Returns `winds`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:179-179`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The buffered wind belongs to the last density callback firing, which for a spacecraft outside the atmosphere may be arbitrarily old or never written at all, in which case zero is reported. The components carry the east, north and vertical convention of the atmosphere driver rather than the planet-frame basis used inside the heat-rate computation, so saved wind is not directly comparable to the relative-velocity terms without applying the same NED rotation.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 63.
