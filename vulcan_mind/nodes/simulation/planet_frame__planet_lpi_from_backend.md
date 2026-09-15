---
id: simulation.planet_frame__planet_lpi_from_backend
label: _planet_lpi_from_backend
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: _planet_lpi_from_backend
  lines:
  - 1
  - 1
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
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: counter
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: true
  description: Positional argument `counter`.
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
  description: Return value of `_planet_lpi_from_backend`.
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

# _planet_lpi_from_backend

## Purpose
`_planet_lpi_from_backend(planet, ephemerides_model, et, counter)` is the uncached path for obtaining the planet-fixed-from-inertial rotation matrix at a given ephemeris time. It exists so that every genuine backend evaluation — typically a SPICE `pxform` frame transformation, which is comparatively expensive — is funnelled through one place where it can be counted for runtime accounting.

## Design & Implementation
Before delegating, it tests `ephemerides_requires_spice(ephemerides_model)` and, when true, performs `Base.Threads.atomic_add!(counter, 1)` on the `Base.Threads.Atomic{Int64}` passed in as `counter`, which the caller supplies as `p.shared_buffers.spice_runtime_counters.planet_pxform_runtime_calls`. It then returns `planet_frame_lpi(planet, et, ephemerides_model)`, declared `::SMatrix{3, 3, Float64}` so the 3x3 direction-cosine matrix stays stack-allocated. The function is `@inline`, so in the common analytic-ephemeris case the predicate folds away and no counting code remains.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `ephemerides_model` | Any | n/a | yes | Positional argument `ephemerides_model`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `counter` | Base.Threads.Atomic{Int64} | n/a | yes | Positional argument `counter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_planet_lpi_from_backend`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl`
- [[simulation.planet_frame__planet_lpi_at|_planet_lpi_at]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:48-48`

**Downstream**

- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:5-5`
<!-- vulcan:connections:end -->

## Limitations
The counter is incremented before the call, so a `planet_frame_lpi` that throws still leaves the count raised — the statistic measures attempts, not successes. Only SPICE-backed models are counted; analytic backends are invisible in the runtime profile even though they are not free. Nothing here validates that the returned matrix is orthonormal or that `et` lies within the loaded kernel coverage; an out-of-coverage time surfaces as a SPICE error from deep inside the backend rather than as a diagnosable message at this boundary.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl` line 1.
