---
id: core.effector_sampling_planetframesample
label: PlanetFrameSample
kind: struct
source:
  file: src/core/types/effector_sampling.jl
  symbol: PlanetFrameSample
  lines:
  - 62
  - 62
inputs:
- id: l_pi
  type: SMatrix{3, 3, Float64, 9}
  units: n/a
  required: true
  description: Field `l_pi`.
- id: pos_pp
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `pos_pp`.
- id: vel_pp
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `vel_pp`.
- id: alt_m
  type: Float64
  units: n/a
  required: true
  description: Field `alt_m`.
- id: lat_rad
  type: Float64
  units: n/a
  required: true
  description: Field `lat_rad`.
- id: lon_rad
  type: Float64
  units: n/a
  required: true
  description: Field `lon_rad`.
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
  type: PlanetFrameSample
  units: n/a
  description: Constructed `PlanetFrameSample`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# PlanetFrameSample

## Purpose
Stage-consistent planet-relative kinematics computed once per ODE evaluation and shared by every effector that requests `planet_frame = true`. It avoids each effector redoing the SPICE frame rotation and geodetic conversion.

## Design & Implementation
Plain immutable struct with `l_pi::SMatrix{3,3,Float64,9}` (rotation from inertial J2000 to planet-fixed IAU frame), `pos_pp` and `vel_pp::SVector{3,Float64}` (planet-relative position in m and velocity in m/s, including the planet rotation term), and scalars `alt_m`, `lat_rad`, `lon_rad`. All fields are computed at the same integrator stage time, which is what "stage-consistent" denotes in the docstring. The struct has no constructor logic and is built by the sampling layer in the simulation engine.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `l_pi` | SMatrix{3, 3, Float64, 9} | n/a | yes | Field `l_pi`. |
| in | `pos_pp` | SVector{3, Float64} | n/a | yes | Field `pos_pp`. |
| in | `vel_pp` | SVector{3, Float64} | n/a | yes | Field `vel_pp`. |
| in | `alt_m` | Float64 | n/a | yes | Field `alt_m`. |
| in | `lat_rad` | Float64 | n/a | yes | Field `lat_rad`. |
| in | `lon_rad` | Float64 | n/a | yes | Field `lon_rad`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlanetFrameSample | n/a | — | Constructed `PlanetFrameSample`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`
- [[simulation.effector_sampling_sample_buffered_planet_frame|sample_buffered_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:152-152`
- [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:48-48`
- [[simulation.effector_sampling_sample_planet_frame_with_lpi|sample_planet_frame_with_lpi]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:57-57`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Whether `alt_m` and `lat_rad` are geodetic or geocentric is not encoded in the type and depends on the engine's `rtolatlong` implementation; effectors that assume one convention can be subtly wrong. `l_pi` is stored as 9 `Float64`s, so the sample is 15 floats plus padding and is copied by value on each hook call. There is no epoch field, so an effector cannot verify that the sample matches the `t` it receives.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 62.
