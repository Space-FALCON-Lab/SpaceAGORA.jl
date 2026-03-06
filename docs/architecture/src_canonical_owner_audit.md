# SpaceAGORA src Canonical-Owner Audit Checklist

Generated on: 2026-03-05T23:28:47.577

## Summary
- Total Julia files under `src/`: 126
- PASS: 126
- FAIL: 0
- Canonical Owner: 114
- Canonical Aggregator: 4
- Scaffold Interface: 8

## File-by-File Checklist
| File | Category | Status | LOC | Wrapper Marker | Silent Stub | Legacy Include Chain |
|---|---|---|---:|---:|---:|---:|
| `src/SpaceAGORA.jl` | Canonical Owner | **PASS** | 30 | no | no | no |
| `src/analysis/plotting/plot_data.jl` | Canonical Owner | **PASS** | 622 | no | no | no |
| `src/analysis/plotting/telemetry_orbit_accuracy_plots.jl` | Canonical Owner | **PASS** | 228 | no | no | no |
| `src/analysis/verification/telemetry_verification.jl` | Canonical Owner | **PASS** | 2335 | no | no | no |
| `src/core/interfaces/reference_system.jl` | Canonical Owner | **PASS** | 428 | no | no | no |
| `src/core/numerics/quaternion_utils.jl` | Canonical Owner | **PASS** | 282 | no | no | no |
| `src/core/simulation_model.jl` | Canonical Aggregator | **PASS** | 92 | no | no | no |
| `src/core/state/reference_system_config.jl` | Canonical Owner | **PASS** | 40 | no | no | no |
| `src/core/state/simulation_configuration.jl` | Canonical Owner | **PASS** | 191 | no | no | no |
| `src/core/types/abstract_types.jl` | Canonical Owner | **PASS** | 11 | no | no | no |
| `src/core/types/runtime_types.jl` | Canonical Owner | **PASS** | 599 | no | no | no |
| `src/core/utils/misc.jl` | Canonical Owner | **PASS** | 21 | no | no | no |
| `src/dynamics/coupled/perturbations.jl` | Canonical Owner | **PASS** | 868 | no | no | no |
| `src/dynamics/models/dynamic_effectors.jl` | Canonical Owner | **PASS** | 28 | no | no | no |
| `src/dynamics/models/thruster_models.jl` | Canonical Owner | **PASS** | 22 | no | no | no |
| `src/dynamics/translational/aerodynamic_models.jl` | Canonical Owner | **PASS** | 966 | no | no | no |
| `src/environment/gravity/gravity_models.jl` | Canonical Owner | **PASS** | 103 | no | no | no |
| `src/environment/atmosphere/density_models.jl` | Canonical Owner | **PASS** | 432 | no | no | no |
| `src/environment/ephemerides/planet_data.jl` | Canonical Owner | **PASS** | 295 | no | no | no |
| `src/environment/ephemerides/planet_shapes.jl` | Canonical Owner | **PASS** | 181 | no | no | no |
| `src/environment/ephemerides/planets.jl` | Canonical Owner | **PASS** | 252 | no | no | no |
| `src/environment/physical_models.jl` | Canonical Owner | **PASS** | 13 | no | no | no |
| `src/vehicle/thermal/thermal_models.jl` | Canonical Owner | **PASS** | 122 | no | no | no |
| `src/vehicle/thermal/thermal_models_module.jl` | Canonical Aggregator | **PASS** | 7 | no | no | no |
| `src/examples/AGORA_Earth.jl` | Canonical Owner | **PASS** | 89 | no | no | no |
| `src/examples/AGORA_Earth_4Sat_Giovanni.jl` | Canonical Owner | **PASS** | 245 | no | no | no |
| `src/examples/AGORA_Earth_Control_Test.jl` | Canonical Owner | **PASS** | 47 | no | no | no |
| `src/examples/AGORA_Earth_GG_Test.jl` | Canonical Owner | **PASS** | 67 | no | no | no |
| `src/examples/AGORA_Earth_SRP_Test.jl` | Canonical Owner | **PASS** | 51 | no | no | no |
| `src/examples/AGORA_Earth_const_torque.jl` | Canonical Owner | **PASS** | 59 | no | no | no |
| `src/examples/AGORA_Keplerian.jl` | Canonical Owner | **PASS** | 89 | no | no | no |
| `src/examples/AGORA_LOFTID.jl` | Canonical Owner | **PASS** | 62 | no | no | no |
| `src/examples/AGORA_Magellan.jl` | Canonical Owner | **PASS** | 62 | no | no | no |
| `src/examples/AGORA_Odyssey.jl` | Canonical Owner | **PASS** | 126 | no | no | no |
| `src/examples/AGORA_Odyssey_Control_Test.jl` | Canonical Owner | **PASS** | 89 | no | no | no |
| `src/examples/AGORA_Titan.jl` | Canonical Owner | **PASS** | 46 | no | no | no |
| `src/examples/AGORA_Titan_Control_Test.jl` | Canonical Owner | **PASS** | 47 | no | no | no |
| `src/examples/AGORA_Vex.jl` | Canonical Owner | **PASS** | 63 | no | no | no |
| `src/examples/AGORA_Vex_Control_Test.jl` | Canonical Owner | **PASS** | 61 | no | no | no |
| `src/examples/CYGNSS_test.jl` | Canonical Owner | **PASS** | 80 | no | no | no |
| `src/examples/Earth_RW_Test.jl` | Canonical Owner | **PASS** | 96 | no | no | no |
| `src/examples/Earth_Thruster_Test.jl` | Canonical Owner | **PASS** | 98 | no | no | no |
| `src/examples/Earth_Torque_Free_Test.jl` | Canonical Owner | **PASS** | 47 | no | no | no |
| `src/examples/GRIFEX_test.jl` | Canonical Owner | **PASS** | 63 | no | no | no |
| `src/examples/typed_example_utils.jl` | Canonical Owner | **PASS** | 151 | no | no | no |
| `src/gnc/control/attitude_control_models.jl` | Canonical Owner | **PASS** | 85 | no | no | no |
| `src/gnc/control/closed_form_solution.jl` | Canonical Owner | **PASS** | 315 | no | no | no |
| `src/gnc/control/control.jl` | Canonical Owner | **PASS** | 295 | no | no | no |
| `src/gnc/control/effectors.jl` | Canonical Owner | **PASS** | 13 | no | no | no |
| `src/gnc/control/eom_ctrl.jl` | Canonical Owner | **PASS** | 351 | no | no | no |
| `src/gnc/control/eoms.jl` | Canonical Owner | **PASS** | 398 | no | no | no |
| `src/gnc/control/heatload_control/second_tsw_calcs.jl` | Canonical Owner | **PASS** | 133 | no | no | no |
| `src/gnc/control/heatload_control/security_mode.jl` | Canonical Owner | **PASS** | 22 | no | no | no |
| `src/gnc/control/heatload_control/time_switch_calcs.jl` | Canonical Owner | **PASS** | 87 | no | no | no |
| `src/gnc/control/heatload_control/utils_timeswitch.jl` | Canonical Owner | **PASS** | 84 | no | no | no |
| `src/gnc/control/legacy_include_helpers.jl` | Canonical Owner | **PASS** | 58 | no | no | no |
| `src/gnc/control/legacy_state_helpers.jl` | Canonical Owner | **PASS** | 40 | no | no | no |
| `src/gnc/control/propulsive_maneuvers.jl` | Canonical Owner | **PASS** | 299 | no | no | no |
| `src/gnc/estimation/estimation.jl` | Scaffold Interface | **PASS** | 16 | no | no | no |
| `src/gnc/guidance/effectors.jl` | Canonical Owner | **PASS** | 12 | no | no | no |
| `src/gnc/guidance/targeting_control/eom_targeting.jl` | Canonical Owner | **PASS** | 854 | no | no | no |
| `src/gnc/guidance/targeting_control/sim_targeting.jl` | Canonical Owner | **PASS** | 523 | no | no | no |
| `src/gnc/guidance/targeting_control/targeting.jl` | Canonical Owner | **PASS** | 292 | no | no | no |
| `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl` | Canonical Owner | **PASS** | 18 | no | no | no |
| `src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl` | Canonical Owner | **PASS** | 4 | no | no | no |
| `src/gnc/navigation/effectors.jl` | Canonical Owner | **PASS** | 11 | no | no | no |
| `src/io/outputs/save_csv.jl` | Canonical Owner | **PASS** | 170 | no | no | no |
| `src/io/outputs/save_results.jl` | Canonical Owner | **PASS** | 408 | no | no | no |
| `src/mission/campaigns/define_mission.jl` | Canonical Owner | **PASS** | 111 | no | no | no |
| `src/mission/campaigns/initial_cond_calc.jl` | Canonical Owner | **PASS** | 32 | no | no | no |
| `src/mission/campaigns/mission_model.jl` | Canonical Owner | **PASS** | 128 | no | no | no |
| `src/mission/campaigns/montecarlo_perturbations.jl` | Canonical Owner | **PASS** | 73 | no | no | no |
| `src/mission/campaigns/montecarlo_set.jl` | Canonical Owner | **PASS** | 86 | no | no | no |
| `src/mission/campaigns/run.jl` | Canonical Owner | **PASS** | 14 | no | no | no |
| `src/mission/campaigns/set_and_run.jl` | Canonical Owner | **PASS** | 24 | no | no | no |
| `src/mission/constellation/network/network.jl` | Scaffold Interface | **PASS** | 20 | no | no | no |
| `src/mission/operations/attitude_control_plans.jl` | Canonical Owner | **PASS** | 340 | no | no | no |
| `src/mission/operations/maneuver_plans.jl` | Canonical Owner | **PASS** | 312 | no | no | no |
| `src/parallel/policy/parallel_policy.jl` | Canonical Owner | **PASS** | 1482 | no | no | no |
| `src/parallel/routing/env_mapping.jl` | Canonical Owner | **PASS** | 183 | no | no | no |
| `src/parallel/routing/outer_route_metrics.jl` | Canonical Owner | **PASS** | 49 | no | no | no |
| `src/parallel/routing/outer_route_selection.jl` | Canonical Owner | **PASS** | 538 | no | no | no |
| `src/parallel/routing/outer_route_state.jl` | Canonical Owner | **PASS** | 198 | no | no | no |
| `src/parallel/routing/parallel_profiles.jl` | Canonical Aggregator | **PASS** | 16 | no | no | no |
| `src/parallel/routing/profile_definitions.jl` | Canonical Owner | **PASS** | 183 | no | no | no |
| `src/simulation/callbacks/callbacks.jl` | Canonical Aggregator | **PASS** | 11 | no | no | no |
| `src/simulation/callbacks/control_callbacks.jl` | Canonical Owner | **PASS** | 55 | no | no | no |
| `src/simulation/callbacks/density_callbacks.jl` | Canonical Owner | **PASS** | 778 | no | no | no |
| `src/simulation/callbacks/event_callbacks.jl` | Canonical Owner | **PASS** | 239 | no | no | no |
| `src/simulation/callbacks/gram_track_cache.jl` | Canonical Owner | **PASS** | 786 | no | no | no |
| `src/simulation/callbacks/navigation_guidance_callbacks.jl` | Canonical Owner | **PASS** | 47 | no | no | no |
| `src/simulation/callbacks/registry.jl` | Canonical Owner | **PASS** | 83 | no | no | no |
| `src/simulation/callbacks/save_fields.jl` | Canonical Owner | **PASS** | 106 | no | no | no |
| `src/simulation/callbacks/thermal_callbacks.jl` | Canonical Owner | **PASS** | 70 | no | no | no |
| `src/simulation/engine/adapters/from_env.jl` | Canonical Owner | **PASS** | 141 | no | no | no |
| `src/simulation/engine/adapters/from_simulation_configuration.jl` | Canonical Owner | **PASS** | 17 | no | no | no |
| `src/simulation/engine/config/artifact_config.jl` | Canonical Owner | **PASS** | 4 | no | no | no |
| `src/simulation/engine/config/parallel_config.jl` | Canonical Owner | **PASS** | 10 | no | no | no |
| `src/simulation/engine/config/runtime_policy_config.jl` | Canonical Owner | **PASS** | 9 | no | no | no |
| `src/simulation/engine/config/simulation_engine_config.jl` | Canonical Owner | **PASS** | 7 | no | no | no |
| `src/simulation/engine/config/solver_config.jl` | Canonical Owner | **PASS** | 9 | no | no | no |
| `src/simulation/engine/dynamics_rhs.jl` | Canonical Owner | **PASS** | 324 | no | no | no |
| `src/simulation/engine/execution.jl` | Canonical Owner | **PASS** | 210 | no | no | no |
| `src/simulation/engine/persistence.jl` | Canonical Owner | **PASS** | 174 | no | no | no |
| `src/simulation/engine/public_api.jl` | Canonical Owner | **PASS** | 5 | no | no | no |
| `src/simulation/engine/reporting.jl` | Canonical Owner | **PASS** | 60 | no | no | no |
| `src/simulation/engine/resume_checkpoint.jl` | Canonical Owner | **PASS** | 56 | no | no | no |
| `src/simulation/engine/setup.jl` | Canonical Owner | **PASS** | 697 | no | no | no |
| `src/simulation/engine/simulation_engine.jl` | Canonical Aggregator | **PASS** | 33 | no | no | no |
| `src/simulation/engine/solver_policy.jl` | Canonical Owner | **PASS** | 354 | no | no | no |
| `src/simulation/events/events.jl` | Canonical Owner | **PASS** | 12 | no | no | no |
| `src/simulation/solver_orchestration/implicit_midpoint_jacobian.jl` | Canonical Owner | **PASS** | 408 | no | no | no |
| `src/simulation/solver_orchestration/include_helpers.jl` | Canonical Owner | **PASS** | 26 | no | no | no |
| `src/simulation/solver_orchestration/integrators.jl` | Canonical Owner | **PASS** | 137 | no | no | no |
| `src/vehicle/actuators/laser_terminal/laser_terminal.jl` | Scaffold Interface | **PASS** | 16 | no | no | no |
| `src/vehicle/actuators/thruster_effectors.jl` | Canonical Owner | **PASS** | 98 | no | no | no |
| `src/vehicle/kinematics/kinematics.jl` | Canonical Owner | **PASS** | 58 | no | no | no |
| `src/vehicle/resources/battery/battery_model.jl` | Scaffold Interface | **PASS** | 16 | no | no | no |
| `src/vehicle/resources/loads/load_model.jl` | Scaffold Interface | **PASS** | 15 | no | no | no |
| `src/vehicle/resources/power_bus/power_bus_model.jl` | Scaffold Interface | **PASS** | 15 | no | no | no |
| `src/vehicle/resources/resources.jl` | Scaffold Interface | **PASS** | 39 | no | no | no |
| `src/vehicle/resources/solar_array/solar_array_model.jl` | Scaffold Interface | **PASS** | 15 | no | no | no |
| `src/vehicle/spacecraft/assembly.jl` | Canonical Owner | **PASS** | 79 | no | no | no |
| `src/vehicle/spacecraft/components.jl` | Canonical Owner | **PASS** | 68 | no | no | no |
| `src/vehicle/spacecraft/model.jl` | Canonical Owner | **PASS** | 259 | no | no | no |
| `src/vehicle/spacecraft/spacecraft_analysis.jl` | Canonical Owner | **PASS** | 335 | no | no | no |

## Failures
- None

## Scaffold Interface Checklist
- [x] `src/gnc/estimation/estimation.jl` (explicit not-implemented errors: yes)
- [x] `src/mission/constellation/network/network.jl` (explicit not-implemented errors: yes)
- [x] `src/vehicle/actuators/laser_terminal/laser_terminal.jl` (explicit not-implemented errors: yes)
- [x] `src/vehicle/resources/battery/battery_model.jl` (explicit not-implemented errors: yes)
- [x] `src/vehicle/resources/loads/load_model.jl` (explicit not-implemented errors: yes)
- [x] `src/vehicle/resources/power_bus/power_bus_model.jl` (explicit not-implemented errors: yes)
- [x] `src/vehicle/resources/resources.jl` (explicit not-implemented errors: yes)
- [x] `src/vehicle/resources/solar_array/solar_array_model.jl` (explicit not-implemented errors: yes)
