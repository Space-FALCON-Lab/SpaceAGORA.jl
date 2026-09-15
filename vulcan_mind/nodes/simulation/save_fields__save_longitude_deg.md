---
id: simulation.save_fields__save_longitude_deg
label: _save_longitude_deg
kind: function
source:
  file: src/simulation/callbacks/save_fields.jl
  symbol: _save_longitude_deg
  lines:
  - 113
  - 113
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
  description: Return value of `_save_longitude_deg`. Returns `longitudes_deg`.
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

# _save_longitude_deg

## Purpose
Save-time getter for planet-fixed longitude in degrees for each spacecraft, supplying the east-west half of the saved ground track.

## Design & Implementation
Marked `@inline` and built like the latitude getter: the ephemeris time `et_start[] + Float64(t)` is computed once, then for each spacecraft `r_intor_p!` rotates the inertial state into the planet-fixed frame and the third element of `rtolatlong` is taken and passed through `rad2deg`. Because the frame is planet-fixed, the saved longitude already includes planetary rotation, so a spacecraft in a fixed inertial orbit shows a steadily drifting longitude in the output.

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
| out | `result` | Any | n/a | — | Return value of `_save_longitude_deg`. Returns `longitudes_deg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:177-177`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:115-115`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:121-121`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:122-122`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:119-119`
<!-- vulcan:connections:end -->

## Limitations
The wrap convention is inherited from `rtolatlong` and is not normalised here, so a consumer plotting a continuous ground track must handle the discontinuity at the branch cut itself. Longitude is degenerate directly over the poles and no guard detects that case. As with the other frame-dependent getters, the transformation is recomputed rather than shared with the altitude and latitude fields.

## Provenance
Mapped from `src/simulation/callbacks/save_fields.jl` line 113.
