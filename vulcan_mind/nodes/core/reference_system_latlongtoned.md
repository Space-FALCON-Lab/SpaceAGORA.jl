---
id: core.reference_system_latlongtoned
label: latlongtoNED
kind: function
source:
  file: src/core/interfaces/reference_system.jl
  symbol: latlongtoNED
  lines:
  - 410
  - 410
inputs:
- id: H_LAN_LON
  type: Any
  units: n/a
  required: true
  description: Positional argument `H_LAN_LON`.
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
  description: Return value of `latlongtoNED`. Returns `SVector{3, SVector{3, Float64}}(uD,
    uN, uE)`.
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

# latlongtoNED

## Purpose
Produces the North, East and Down unit vectors of the local frame at a latitude and longitude, expressed in planet-fixed axes.

## Design & Implementation
Forms the three unit vectors in an intermediate frame whose x-z plane contains the position, then rotates them about z by the longitude with a static rotation matrix. Returns an `SVector` of three `SVector`s in the order down, north, east. The input is read as `[h, lat, lon]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `H_LAN_LON` | Any | n/a | yes | Positional argument `H_LAN_LON`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `latlongtoNED`. Returns `SVector{3, SVector{3, Float64}}(uD, uN, uE)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:748-748`
- [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:391-391`
- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:126-126`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:138-138`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:155-155`
- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:73-73`
- [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:344-344`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:160-160`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:126-126`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:138-138`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:155-155`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:160-160`
- [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_out` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:748-748`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:126-126`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:155-155`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:37-37`
- [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:37-37`
- [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:242-242`

**Downstream**

- `callees` → [[core.reference_system_orbital_elements_to_lvlh_quaternion|orbital_elements_to_lvlh_quaternion]] · `callers` · call · `src/core/interfaces/reference_system.jl:449-449`
<!-- vulcan:connections:end -->

## Limitations
The return order is down, north, east while the function name says NED, so a caller destructuring positionally as N, E, D gets the wrong axes; the input takes altitude first although altitude is unused.

## Provenance
Mapped from `src/core/interfaces/reference_system.jl` line 410.
