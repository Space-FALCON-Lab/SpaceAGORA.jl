---
id: core.reference_system_rtolatlong
label: rtolatlong
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: rtolatlong
  lines:
  - 351
  - 351
inputs:
- id: r_p
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r_p`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: spherical_harmonic_topography
  type: Bool
  units: n/a
  required: false
  description: Positional argument `spherical_harmonic_topography` (default `false`).
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
  type: SVector
  units: n/a
  description: Return value of `rtolatlong`. Returns `SVector{3, Float64}([alt, lat,
    lon])`.
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

# rtolatlong

## Purpose
Converts a planet-fixed position to geodetic altitude, latitude and longitude by Bowring's closed-form method, optionally measuring altitude against a spherical-harmonic topography.

## Theory & Math
$$
\theta = \operatorname{atan2}(z R_e,\; p R_p),\qquad \phi = \operatorname{atan2}\left(z + e'^2 R_p \sin^3\theta,\; p - e^2 R_e \cos^3\theta\right)
$$

with $p = \sqrt{x^2 + y^2}$, $e^2 = 1 - (1-f)^2$ and $e'^2 = e^2 / (1 - e^2)$.

## Design & Implementation
Computes flattening, first and second eccentricity squared, and the equatorial distance. Bowring's auxiliary angle `θ = atan(z R_e, p R_p)` gives geodetic latitude in one step without iteration; longitude is `atan(y, x)`. Altitude defaults to Bowring's closed form using the prime-vertical radius, or, with the topography flag, to the radial distance minus a `planet.topography_function` evaluated with the including module's `args` global if defined. A three-argument method accepting an ephemerides model simply forwards. Returns `[alt, lat, lon]` as an `SVector`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_p` | SVector{3, Float64} | n/a | yes | Positional argument `r_p`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `spherical_harmonic_topography` | Bool | n/a | no | Positional argument `spherical_harmonic_topography` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `rtolatlong`. Returns `SVector{3, Float64}([alt, lat, lon])`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__telemetry_altitude_km|_telemetry_altitude_km]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:71-71`
- [[analysis.telemetry_loading__extract_extrema_series|_extract_extrema_series]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:96-96`
- [[dynamics.aerodynamic_wrench_models__aero_link_atmosphere_query|_aero_link_atmosphere_query]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:507-507`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2055-2055`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1910-1910`
- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:127-127`
- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:120-120`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:132-132`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:149-149`
- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:71-71`
- [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:334-334`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:129-129`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:120-120`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:132-132`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:149-149`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:129-129`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:120-120`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:127-127`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/save_fields.jl:94-94`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/effector_sampling.jl:47-47`
- [[simulation.effector_sampling_sample_planet_frame|sample_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:47-47`
- [[simulation.effector_sampling_sample_planet_frame_with_lpi|sample_planet_frame_with_lpi]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:56-56`
- [[simulation.save_fields__save_altitude|_save_altitude]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:94-94`
- [[simulation.save_fields__save_latitude_deg|_save_latitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:108-108`
- [[simulation.save_fields__save_longitude_deg|_save_longitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:122-122`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:105-105`
- [[simulation.targeting__gram_linear_target|_gram_linear_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:181-181`
- [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:272-272`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:237-237`
- [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:241-241`
- [[simulation_a.vacuum_predicted_gram_build_vacuum_gram_cache__build_vacuum_gram_cache_bang|_build_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:213-213`

**Downstream**

- `callees` → [[core.reference_system__planet_flattening|_planet_flattening]] · `callers` · call · `src/core/interfaces/reference_system.jl:357-357`
<!-- vulcan:connections:end -->

## Limitations
The topography path reaches for a module-level `args` global by name, a duck-typed contract that only the sandbox test exercises; the returned triple orders altitude first, which differs from the input order of `latlongtor`.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 351.
