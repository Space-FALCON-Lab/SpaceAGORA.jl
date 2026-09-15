---
id: core.reference_system_r_intor_p_bang
label: r_intor_p!
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: r_intor_p!
  lines:
  - 26
  - 26
inputs:
- id: r_i
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r_i`.
- id: v_i
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `v_i`.
- id: planet
  type: T
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: Tuple{SVector{3,
  units: n/a
  description: 'Return value of `r_intor_p!`; mutates `r_i` in place. Type parameters:
    `T`.'
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

# r_intor_p!

## Purpose
Rotates an inertial position and velocity into the planet-centred, planet-fixed frame, accounting for the frame's rotation in the velocity.

## Theory & Math
$$
r_p = L_{PI}\, r_i,\qquad v_p = L_{PI}\, v_i - \omega \times r_p
$$

## Design & Implementation
Three methods. The legacy two-vector form uses the planet's stored `L_PI` matrix and subtracts the transport term `ω × r_p` after rotating, because `planet.ω` is expressed in planet-fixed axes. The `et` form delegates to `_j2000_to_body_fixed_state` for a SPICE state transform. The `ephemerides_model` form chooses between the SPICE route and an analytic `planet_frame_lpi` rotation with the same transport-term subtraction. Despite the `!` all methods are pure and return a fresh tuple of static vectors.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `r_i` | SVector{3, Float64} | n/a | yes | Positional argument `r_i`. |
| in | `v_i` | SVector{3, Float64} | n/a | yes | Positional argument `v_i`. |
| in | `planet` | T | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `r_intor_p!`; mutates `r_i` in place. Type parameters: `T`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__telemetry_altitude_km|_telemetry_altitude_km]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:70-70`
- [[analysis.telemetry_loading__extract_extrema_series|_extract_extrema_series]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:95-95`
- [[core.reference_system__planet_flattening|_planet_flattening]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:91-91`
- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:603-603`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:912-912`
- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:125-125`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:99-99`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:116-116`
- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:70-70`
- [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:333-333`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:81-81`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:44-44`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:45-45`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:81-81`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:99-99`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:125-125`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/save_fields.jl:93-93`
- [[simulation.save_fields__save_altitude|_save_altitude]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:93-93`
- [[simulation.save_fields__save_latitude_deg|_save_latitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:107-107`
- [[simulation.save_fields__save_longitude_deg|_save_longitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:121-121`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:104-104`
- [[simulation.targeting__gram_linear_target|_gram_linear_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:180-180`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:236-236`
- [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:240-240`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The name promises mutation that does not occur, which is misleading; the legacy form uses a fixed `L_PI` so it ignores the planet's rotation since epoch and is only correct at the epoch the matrix was computed for.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 26.
