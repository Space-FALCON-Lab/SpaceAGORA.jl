---
id: simulation.save_fields__save_velocities
label: _save_velocities
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_velocities
  lines:
  - 25
  - 25
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
  description: Return value of `_save_velocities`. Returns `velocities`.
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

# _save_velocities

## Purpose
Save-time getter for spacecraft inertial velocities, returning a `Vector{SVector{3, Float64}}` of length `num_sats` in metres per second, one entry per spacecraft.

## Design & Implementation
Structurally the mirror of the position getter: `@inline`, an `undef` result vector of static 3-vectors, and an `@inbounds` loop over `_simulation_engine_module()._state_velocity_ii(u, i)`. Using `SVector{3, Float64}` rather than a plain `Vector` keeps each sample stack-allocated and avoids a second level of heap indirection per spacecraft. The signature carries `t` and `integrator` only for getter uniformity; neither is read, since velocity is a primary state component and needs no time-dependent transformation.

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
| out | `result` | Any | n/a | — | Return value of `_save_velocities`. Returns `velocities`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:174-174`

**Downstream**

- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:28-28`
<!-- vulcan:connections:end -->

## Limitations
Like the other state getters it allocates one vector per save point. It reports the raw inertial velocity, not the planet-relative or wind-relative velocity that the aerodynamic computations use, so a consumer comparing saved velocity against drag must apply the frame and wind corrections itself. No bounds check protects against `num_sats` exceeding what the state holds.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 25.
