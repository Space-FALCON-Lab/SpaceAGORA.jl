---
id: environment.simple_ephemerides_planet_frame_lpi
label: planet_frame_lpi
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: planet_frame_lpi
  lines:
  - 92
  - 92
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
  description: Return value of `planet_frame_lpi`.
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

# planet_frame_lpi

## Purpose
Provides the J2000-to-planet-centred-planet-fixed rotation matrix at time `et` for either ephemeris model, so that geographic models (atmosphere, gravity harmonics) can be evaluated in body-fixed coordinates.

## Design & Implementation
The `SpiceEphemeridesModel` method delegates to `_spice_planet_frame_lpi(planet, et)`. The `SimpleEphemeridesModel` method computes `elapsed = et - model.reference_epoch_seconds` and the spin angle `θ`: when `prime_meridian_at_reference_rad` is `NaN` (planet-true default) Earth uses `_earth_gmst_iau82_rad(elapsed)` and other planets use `planet.ω[3] * elapsed`; an explicit finite value selects the legacy linear model `pm + planet.ω[3] * elapsed`. The angle is passed to `_rotation_about_spin_axis(θ)` and an `SMatrix{3,3,Float64}` is returned. Both methods are `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SMatrix{3, | n/a | — | Return value of `planet_frame_lpi`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:511-511`
- [[core.reference_system__planet_flattening|_planet_flattening]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:108-108`
- [[dynamics.perturbations__harmonics_lpi_at_bang|_harmonics_lpi_at!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:588-588`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2052-2052`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1907-1907`
- [[environment.simple_ephemerides__earth_gmst_iau82_rad|_earth_gmst_iau82_rad]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:105-105`
- [[gnc.targeting_control__edg_planet_frame_lpi|_edg_planet_frame_lpi]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:325-325`
- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:180-180`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- [[simulation.planet_frame__planet_lpi_from_backend|_planet_lpi_from_backend]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:5-5`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1885-1885`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:221-221`
- [[vehicle.model__initial_condition_lpi|_initial_condition_lpi]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:99-99`

**Downstream**

- `callees` → [[environment.simple_ephemerides__spice_planet_frame_lpi|_spice_planet_frame_lpi]] · `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:93-93`
<!-- vulcan:connections:end -->

## Limitations
The simple model treats the spin axis as J2000 z and the rate `ω[3]` as constant, ignoring precession, nutation, and polar motion. Non-Earth planets in default mode start at zero prime-meridian angle at the reference epoch, which is not their true IAU orientation, so absolute longitudes are only meaningful relative to that epoch. The `NaN` sentinel means an accidental NaN from upstream arithmetic silently switches modes.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 92.
