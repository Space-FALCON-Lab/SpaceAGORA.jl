---
id: environment.simple_ephemerides__spice_planet_frame_lpi
label: _spice_planet_frame_lpi
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: _spice_planet_frame_lpi
  lines:
  - 79
  - 79
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
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
  description: Return value of `_spice_planet_frame_lpi`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# _spice_planet_frame_lpi

## Purpose
Obtains the J2000-to-planet-fixed rotation matrix at ephemeris time `et` from SPICE for the `SpiceEphemeridesModel`, with an automatic fallback from the high-precision ITRF93 Earth frame to IAU_EARTH when the required kernel is missing.

## Design & Implementation
All SPICE calls are wrapped in `lock(SPICE_LOCK) do ... end` because the CSPICE library is not thread-safe. For `planet.name == "Earth"` it attempts `pxform("J2000", "ITRF93", et)` inside a `try` block and, on any exception, retries with `"IAU_EARTH"`. For all other bodies it calls `pxform("J2000", _spice_body_fixed_frame(planet), et)` once. Each result is wrapped as `SMatrix{3,3,Float64}`. The `return` inside the `do` block returns from the closure and thereby the function.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `_spice_planet_frame_lpi`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:93-93`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`

**Downstream**

- `callees` → [[core.reference_system__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:88-88`
- `callees` → [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:88-88`
- `callees` → [[environment.simple_ephemerides__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:88-88`
<!-- vulcan:connections:end -->

## Limitations
The bare `catch` swallows every exception type, including unrelated errors such as an invalid `et`, and silently degrades Earth fidelity by up to tens of metres at the surface without logging. The `try` incurs a full SPICE error-handling round trip on every call when ITRF93 is absent, which is slow inside an integrator right-hand side. Holding the global lock serialises all ephemeris evaluations across threads.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 79.
