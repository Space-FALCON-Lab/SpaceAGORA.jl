---
id: vehicle.model__initial_condition_lpi
label: _initial_condition_lpi
kind: function
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: _initial_condition_lpi
  lines:
  - 92
  - 92
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: L_PI
  type: Any
  units: n/a
  required: true
  description: Positional argument `L_PI`.
- id: initial_time
  type: Any
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: ephemerides_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `ephemerides_model`.
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `_initial_condition_lpi`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# _initial_condition_lpi

## Purpose
Resolves which inertial-to-planet-fixed rotation `L_PI` the oblate initial-condition constructor should use, with a fixed precedence over an explicit matrix, an ephemeris evaluation at `initial_time`, a matrix stored on the planet, and finally identity.

## Design & Implementation
`@inline _initial_condition_lpi(planet, L_PI, initial_time, ephemerides_model)::SMatrix{3,3,Float64,9}`. If `L_PI !== nothing` it is converted to an `SMatrix` and returned. Else if `initial_time !== nothing`, the ephemeris model (defaulting to `SpiceEphemeridesModel()` when `ephemerides_model === nothing`) is used to compute `planet_frame_lpi(planet, ephemerides_time_seconds(initial_time, model), model)`. Else if `hasproperty(planet, :L_PI)` and that matrix is non-zero (`sum(abs2, ...) > 0`), it is returned. Otherwise the constant identity `I3` is returned, meaning the inertial and planet-fixed frames are assumed aligned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `L_PI` | Any | n/a | yes | Positional argument `L_PI`. |
| in | `initial_time` | Any | n/a | yes | Positional argument `initial_time`. |
| in | `ephemerides_model` | Any | n/a | yes | Positional argument `ephemerides_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_initial_condition_lpi`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`
- [[vehicle.model__initial_condition_radius_for_oblate_altitude|_initial_condition_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:201-201`

**Downstream**

- `callees` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/vehicle/spacecraft/model.jl:99-99`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/vehicle/spacecraft/model.jl:99-99`
- `callees` → [[environment.simple_ephemerides_spiceephemeridesmodel|SpiceEphemeridesModel]] · `callers` · call · `src/vehicle/spacecraft/model.jl:98-98`
<!-- vulcan:connections:end -->

## Limitations
Falling back to identity silently changes the meaning of apsis altitudes on an oblate planet; no warning is issued. A planet `L_PI` of all zeros is treated as uninitialised, which is a heuristic rather than an explicit flag. Constructing a default `SpiceEphemeridesModel()` requires SPICE kernels to be loaded, and errors surface from inside `planet_frame_lpi`.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 92.
