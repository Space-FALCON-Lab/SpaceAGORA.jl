---
id: simulation.save_fields__save_altitude
label: _save_altitude
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_altitude
  lines:
  - 85
  - 85
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
  description: Return value of `_save_altitude`. Returns `altitudes`.
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

# _save_altitude

## Purpose
Save-time getter for geometric altitude above the planet surface, in metres, for each spacecraft at the current save point.

## Design & Implementation
Marked `@inline`. It forms the ephemeris time as `integrator.p.shared_buffers.et_start[] + Float64(t)`, hoisting that computation and the `planet` and `ephemerides_model` lookups out of the per-spacecraft loop. For each spacecraft it reads inertial position and velocity, rotates into the planet-fixed frame with `r_intor_p!(pos, vel, planet, et, ephemerides_model)`, and takes the first element of `rtolatlong(rp, planet, ephemerides_model)`, which is altitude. The second returned value of `r_intor_p!` is discarded with `_`.

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
| out | `result` | Any | n/a | — | Return value of `_save_altitude`. Returns `altitudes`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:175-175`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:87-87`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:93-93`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:94-94`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:91-91`
<!-- vulcan:connections:end -->

## Limitations
The bang in `r_intor_p!` signals that the rotation routine mutates shared workspace, so this getter is not safe to run concurrently with another consumer of that workspace. The per-spacecraft `rtolatlong` call allocates and returns all three of altitude, latitude and longitude while only altitude is kept, and the sibling latitude and longitude getters repeat the identical frame conversion, so a save point with all three fields enabled performs the transformation three times per spacecraft.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 85.
