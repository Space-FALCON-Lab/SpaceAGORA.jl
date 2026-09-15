---
id: grp.src_gnc_guidance
label: gnc/guidance/
kind: group
inputs:
- id: members_in
  type: call
  units: n/a
  required: false
  description: Calls into any member of this block from outside it.
outputs:
- id: members_out
  type: call
  units: n/a
  description: Calls from any member of this block to nodes outside it.
expands: module.gnc
tags:
- cluster
charts:
- gnc
origin: agent
---

# gnc/guidance/

## Purpose
The guidance layer: decides what the vehicle should do next — aerobraking pass strategy, target energy, RPO paths, propulsive manoeuvres — and hands commands to control.

## Design & Implementation
`aerobraking/` holds the E-EDG and T-EDG strategies behind a dispatcher; `target_energy_bracketing.jl` is the energy-depletion guidance model; `rpo/` is the HYPR path planner and its comparison baselines; `thruster_guidance/` schedules campaign burns; `guidance_hooks.jl` and `guidance_models.jl` define the interface.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `members_in` | call | n/a | no | Calls into any member of this block from outside it. |
| out | `members_out` | call | n/a | — | Calls from any member of this block to nodes outside it. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.assets|RPOStationAssets — internals]] · `members_out` → `members_in` · call · `src/assets/rpo_station_assets.jl:88-88`
- [[grp.cli|SpaceAGORACLI — internals]] · `members_out` → `members_in` · call · `src/cli/assets.jl:52-52`
- [[grp.mission|AerobrakingPolicy — internals]] · `members_out` → `members_in` · call · `src/mission/operations/maneuver_plans.jl:40-40`
- [[grp.src_analysis_verification|analysis/verification/]] · `members_out` → `members_in` · call · `src/analysis/verification/telemetry_verification/calibration.jl:84-84`
- [[grp.src_analysis_visualization|analysis/visualization/]] · `members_out` → `members_in` · call · `src/analysis/visualization/rpo/rpo_visualization.jl:10-10`
- [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_out` → `members_in` · call · `src/dynamics/coupled/perturbations.jl:886-886`
- [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_out` → `members_in` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:303-303`
- [[grp.src_environment_ephemerides|environment/ephemerides/]] · `members_out` → `members_in` · call · `src/environment/ephemerides/planets.jl:47-47`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:462-462`
- [[grp.src_gnc_hypr|gnc/hypr/]] · `members_out` → `members_in` · call · `src/gnc/hypr/hypr_utils.jl:160-160`
- [[grp.src_gnc_robotics|gnc/robotics/]] · `members_out` → `members_in` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:289-289`
- [[grp.src_parallel_policy|parallel/policy/]] · `members_out` → `members_in` · call · `src/parallel/policy/context.jl:302-302`
- [[grp.src_parallel_routing|parallel/routing/]] · `members_out` → `members_in` · call · `src/parallel/routing/outer_route_selection.jl:409-409`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `members_in` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:148-148`
- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `members_in` · call · `src/simulation/campaigns/adaptive_routing.jl:197-197`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `members_in` · call · `src/simulation/engine/execution.jl:78-78`
- [[grp.src_vehicle_robotics|vehicle/robotics/]] · `members_out` → `members_in` · call · `src/vehicle/robotics/robotics.jl:189-189`
- [[grp.src_vehicle_spacecraft|vehicle/spacecraft/]] · `members_out` → `members_in` · call · `src/vehicle/spacecraft/assembly.jl:59-59`
- [[grp.src_vehicle_structure|vehicle/structure/]] · `members_out` → `members_in` · call · `src/vehicle/structure/mass_properties.jl:105-105`
- [[module.gnc|CommandTypes]] · `api` → `members_in` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl`

**Downstream**

- `members_out` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- `members_out` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1041-1041`
- `members_out` → [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:11-11`
- `members_out` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:155-155`
- `members_out` → [[core.reference_system_orbitalelemtorv|orbitalelemtorv]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:117-117`
- `members_out` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:125-125`
- `members_out` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:127-127`
- `members_out` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_constant|aerodynamic_coefficient_constant]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:316-316`
- `members_out` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:195-195`
- `members_out` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_no_ballistic_flight|aerodynamic_coefficient_no_ballistic_flight]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:320-320`
- `members_out` → [[dynx.coupled_perturbations_srp|srp]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:256-256`
- `members_out` → [[environment.density_models_density_polyfit|density_polyfit]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:191-191`
- `members_out` → [[environment.gravity_models_aerobraking_gravity_force_ii|aerobraking_gravity_force_ii]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:236-236`
- `members_out` → [[grp.cli|SpaceAGORACLI — internals]] · `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1041-1041`
- `members_out` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- `members_out` → [[grp.src_core_interfaces|core/interfaces/]] · `members_in` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:117-117`
- `members_out` → [[grp.src_core_state|core/state/]] · `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:91-91`
- `members_out` → [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_in` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:195-195`
- `members_out` → [[grp.src_environment_gravity|environment/gravity/]] · `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:236-236`
- `members_out` → [[grp.src_gnc_command_types_jl|gnc/command_types.jl]] · `members_in` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:122-122`
- `members_out` → [[grp.src_gnc_control|gnc/control/]] · `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:214-214`
- `members_out` → [[grp.src_gnc_hypr|gnc/hypr/]] · `members_in` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:50-50`
- `members_out` → [[grp.src_gnc_internal|gnc/internal/]] · `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1037-1037`
- `members_out` → [[grp.src_gnc_navigation|gnc/navigation/]] · `members_in` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:5-5`
- `members_out` → [[grp.src_gnc_robotics|gnc/robotics/]] · `members_in` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:368-368`
- `members_out` → [[grp.src_vehicle_structure|vehicle/structure/]] · `members_in` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:156-156`
- `members_out` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:161-161`
- `members_out` → [[vehicle.thermal_models_heatrate_convective_radiative|heatrate_convective_radiative]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:228-228`
<!-- vulcan:connections:end -->

## Limitations
Guidance and control for aerobraking are tightly coupled through shared state structs rather than through the command types alone.

## Provenance
Macro block generated by `vulcan compile` (D13) from `src/gnc/guidance`; groups 265 nodes.
Its members are listed in chart `gnc-gnc-guidance`.
