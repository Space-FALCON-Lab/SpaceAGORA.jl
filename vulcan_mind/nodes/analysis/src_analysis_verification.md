---
id: grp.src_analysis_verification
label: analysis/verification/
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
expands: module.analysis
tags:
- cluster
charts:
- analysis
origin: agent
---

# analysis/verification/

## Purpose
The verification half of the analysis module: everything that compares a SpaceAGORA simulation against flight telemetry and scores the match.

## Design & Implementation
Contains the `telemetry_verification.jl` entry file and the `telemetry_verification/` directory beneath it. A manifest of scenarios enters; per-scenario simulation configurations are built, run, compared channel by channel against loaded telemetry, and summarised into CSV tables and plots. It is the only consumer of the engine's results that feeds back a pass/fail verdict.

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

- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `members_in` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:138-138`
- [[ext.spaceagoragramsuiteext__gram_utc_string|_gram_utc_string]] · `callees` → `members_in` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:68-68`
- [[grp.assets|RPOStationAssets — internals]] · `members_out` → `members_in` · call · `src/assets/rpo_station_assets.jl:119-119`
- [[grp.io|IOConfig — internals]] · `members_out` → `members_in` · call · `src/io/serialization/io_serialization.jl:75-75`
- [[grp.mission|AerobrakingPolicy — internals]] · `members_out` → `members_in` · call · `src/mission/operations/maneuver_plans.jl:2-2`
- [[grp.src_core_numerics|core/numerics/]] · `members_out` → `members_in` · call · `src/core/numerics/quaternion_utils.jl:20-20`
- [[grp.src_core_state|core/state/]] · `members_out` → `members_in` · call · `src/core/state/simulation_configuration.jl:223-223`
- [[grp.src_core_types|core/types/]] · `members_out` → `members_in` · call · `src/core/types/effector_sampling.jl:51-51`
- [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_out` → `members_in` · call · `src/dynamics/coupled/perturbations.jl:1518-1518`
- [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_out` → `members_in` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:183-183`
- [[grp.src_dynamics_rotational|dynamics/rotational/]] · `members_out` → `members_in` · call · `src/dynamics/rotational/torque_models.jl:6-6`
- [[grp.src_dynamics_translational|dynamics/translational/]] · `members_out` → `members_in` · call · `src/dynamics/translational/point_mass_dynamics.jl:19-19`
- [[grp.src_environment_atmosphere|environment/atmosphere/]] · `members_out` → `members_in` · call · `src/environment/atmosphere/density_models.jl:840-840`
- [[grp.src_environment_ephemerides|environment/ephemerides/]] · `members_out` → `members_in` · call · `src/environment/ephemerides/simple_ephemerides.jl:49-49`
- [[grp.src_environment_gravity|environment/gravity/]] · `members_out` → `members_in` · call · `src/environment/gravity/gravity_models.jl:87-87`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `members_in` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `members_in` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- [[grp.src_gnc_hypr|gnc/hypr/]] · `members_out` → `members_in` · call · `src/gnc/hypr/hypr_utils.jl:32-32`
- [[grp.src_gnc_internal|gnc/internal/]] · `members_out` → `members_in` · call · `src/gnc/internal/bridge_helpers.jl:110-110`
- [[grp.src_gnc_navigation|gnc/navigation/]] · `members_out` → `members_in` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:18-18`
- [[grp.src_gnc_robotics|gnc/robotics/]] · `members_out` → `members_in` · call · `src/gnc/robotics/robot_arm_hypr/config.jl:8-8`
- [[grp.src_parallel_policy|parallel/policy/]] · `members_out` → `members_in` · call · `src/parallel/policy/persistent_hints.jl:235-235`
- [[grp.src_parallel_routing|parallel/routing/]] · `members_out` → `members_in` · call · `src/parallel/routing/outer_route_selection.jl:472-472`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `members_in` · call · `src/simulation/callbacks/event_callbacks.jl:151-151`
- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `members_in` · call · `src/simulation/campaigns/adaptive_routing.jl:75-75`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `members_in` · call · `src/simulation/engine/dynamics_rhs.jl:178-178`
- [[grp.src_vehicle_robotics|vehicle/robotics/]] · `members_out` → `members_in` · call · `src/vehicle/robotics/robotics.jl:315-315`
- [[grp.src_vehicle_spacecraft|vehicle/spacecraft/]] · `members_out` → `members_in` · call · `src/vehicle/spacecraft/model.jl:198-198`
- [[io.io_serialization__load_checkpoint|_load_checkpoint]] · `callees` → `members_in` · call · `src/io/serialization/io_serialization.jl:75-75`
- [[module.analysis|TelemetryVerification]] · `api` → `members_in` · call · `src/analysis/verification/telemetry_verification/calibration.jl`

**Downstream**

- `members_out` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:184-184`
- `members_out` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:84-84`
- `members_out` → [[grp.cli|SpaceAGORACLI — internals]] · `members_in` · call · `src/analysis/verification/telemetry_verification/runner.jl:361-361`
- `members_out` → [[grp.src_core_interfaces|core/interfaces/]] · `members_in` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:70-70`
- `members_out` → [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_in` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:161-161`
- `members_out` → [[grp.src_gnc_guidance|gnc/guidance/]] · `members_in` · call · `src/analysis/verification/telemetry_verification/calibration.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations
It verifies orbital state channels only; attitude, thermal and power telemetry have no comparison path here.

## Provenance
Macro block generated by `vulcan compile` (D13) from `src/analysis/verification`; groups 159 nodes.
Its members are listed in chart `analysis-analysis-verification`.
