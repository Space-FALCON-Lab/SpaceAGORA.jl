---
id: simulation.planet_frame__planet_lpi_at
label: _planet_lpi_at
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: _planet_lpi_at
  lines:
  - 40
  - 40
inputs:
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
  type: SMatrix{3,
  units: n/a
  description: Return value of `_planet_lpi_at`.
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

# _planet_lpi_at

## Purpose
`_planet_lpi_at(p, t::Float64)` is the single entry point that turns a simulation-relative time into the planet-fixed-from-inertial rotation matrix, choosing transparently between the interpolation cache and the ephemeris backend. Everything that needs the body-fixed frame during a step — the density callback, the aerodynamics, the topography lookup — goes through it.

## Design & Implementation
It pulls `planet` and `ephemerides_model` out of `p.args.environment_model`, converts the relative time to absolute ephemeris time with `et = p.shared_buffers.et_start[] + t`, and picks up the SPICE call counter `p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls`. It then reads `p.shared_buffers.planet_frame_ephemeris_cache[]` from its `Ref`; if the stored value `isa PlanetFrameEphemerisCache` it tries `_planet_lpi_from_cache(cache_entry, et)` and falls back to `_planet_lpi_from_backend` only when that returns `nothing`, otherwise it goes straight to the backend. The `isa` test rather than a `!== nothing` test means an unset or differently-typed cache slot degrades safely. The return type is annotated `::SMatrix{3, 3, Float64}` and the function is `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_planet_lpi_at`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.effector_sampling__planet_lpi_at_engine|_planet_lpi_at_engine]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:39-39`
- [[simulation.planet_frame_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:73-73`
- [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:270-270`
- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:73-73`
- [[simulation_a.vacuum_predicted_gram_build_vacuum_gram_cache__build_vacuum_gram_cache_bang|_build_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:211-211`

**Downstream**

- `callees` → [[simulation.planet_frame__planet_lpi_from_backend|_planet_lpi_from_backend]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:48-48`
- `callees` → [[simulation.planet_frame__planet_lpi_from_cache|_planet_lpi_from_cache]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
`et_start[]` must already have been populated by the callback's initialiser; if `_planet_lpi_at` is called before initialisation, `et` is simply `t` seconds from the J2000 epoch and the resulting frame is silently wrong rather than an error. The cache `Ref` is read without synchronisation, so a concurrent refresh from another task can be observed mid-swap. Every cache miss costs the full bracketing search before the backend call, so a trajectory that repeatedly steps just outside the cached window pays for both paths on every step.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl` line 40.
