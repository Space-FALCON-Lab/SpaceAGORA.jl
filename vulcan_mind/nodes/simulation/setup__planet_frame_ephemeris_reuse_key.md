---
id: simulation.setup__planet_frame_ephemeris_reuse_key
label: _planet_frame_ephemeris_reuse_key
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _planet_frame_ephemeris_reuse_key
  lines:
  - 289
  - 289
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: ephemerides_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `ephemerides_model`.
- id: et_start
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et_start`.
- id: mission_end_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_end_s`.
- id: dt_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `dt_s`.
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
  type: PlanetFrameEphemerisReuseKey
  units: n/a
  description: Return value of `_planet_frame_ephemeris_reuse_key`.
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

# _planet_frame_ephemeris_reuse_key

## Purpose
Produces the reuse key for a `PlanetFrameEphemerisCache`, capturing everything that determines the tabulated rotation: planet identity and constants, ephemerides model, and time grid.

## Design & Implementation
Takes `planet`, `ephemerides_model`, `et_start`, `mission_end_s`, `dt_s` and returns the 6-tuple `(string(planet.name), _ephemerides_model_reuse_key(ephemerides_model), _planet_transform_key(planet), _cache_time_key(et_start), _cache_time_key(mission_end_s), _cache_time_key(dt_s))` of type `PlanetFrameEphemerisReuseKey = Tuple{String, String, NTuple{9, Int64}, Int64, Int64, Int64}`. Including `_planet_transform_key` means editing a planet's radius or pole in a scenario invalidates the cache automatically.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `ephemerides_model` | Any | n/a | yes | Positional argument `ephemerides_model`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `dt_s` | Float64 | n/a | yes | Positional argument `dt_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PlanetFrameEphemerisReuseKey | n/a | — | Return value of `_planet_frame_ephemeris_reuse_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1872-1872`

**Downstream**

- `callees` → [[simulation.setup__cache_time_key|_cache_time_key]] · `callers` · call · `src/simulation/engine/setup.jl:294-294`
- `callees` → [[simulation.setup__ephemerides_model_reuse_key|_ephemerides_model_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:292-292`
- `callees` → [[simulation.setup__planet_transform_key|_planet_transform_key]] · `callers` · call · `src/simulation/engine/setup.jl:293-293`
<!-- vulcan:connections:end -->

## Limitations
`string(planet.name)` accepts a `Symbol` or `String` but two planets with the same name and constants but different, unkeyed properties (for example a custom rotation model) would collide. The μ quantisation in `_planet_transform_key` is coarse. Keys are large tuples hashed on every lookup, which is negligible at setup but not free.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 289.
