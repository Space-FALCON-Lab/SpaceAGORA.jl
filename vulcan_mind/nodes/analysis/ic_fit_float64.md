---
id: analysis.ic_fit_float64
label: Float64
kind: function
source:
  file: src/analysis/verification/telemetry_verification/ic_fit.jl
  symbol: Float64
  lines:
  - 24
  - 24
inputs:
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
  type: Any
  units: n/a
  description: Return value of `Float64`. Returns `> (Float64(sub.sim_interp_value_km[i]),
    Float64(sub.error_km[i]))`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# Float64

## Purpose
This node was extracted from the `Float64(...)` conversions on line 24 of `ic_fit.jl`, inside `_ic_fit_series`. It is not a user-defined function: it records the point where per-sample telemetry CSV columns are coerced to `Float64` while building the axis-keyed residual series used by the differential-correction fit `fit_initial_state`.

## Design & Implementation
`_ic_fit_series(errors_csv, scenario_name)` reads the errors CSV into a `DataFrame`, filters rows to the scenario and to the three events in `_IC_FIT_EVENTS` (`state_x_time`, `state_y_time`, `state_z_time`), and throws an `ArgumentError` if no rows remain. For each event it builds a `Dict{Float64, NTuple{2,Float64}}` mapping `Float64(telemetry_axis)` (the comparison time key) to `(Float64(sim_interp_value_km), Float64(error_km))`. Those dictionaries are later intersected across the baseline and six perturbation runs so the finite-difference Jacobian in `fit_initial_state` uses only common sample times.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `Float64`. Returns `> (Float64(sub.sim_interp_value_km[i]), Float64(sub.error_km[i]))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.calibration__calibration_score|_calibration_score]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:61-61`
- [[analysis.decay_diagnostics_flight_density_table|flight_density_table]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:158-158`
- [[analysis.manifest_parsing__optional_float64_vector|_optional_float64_vector]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:139-139`
- [[analysis.manifest_parsing__optional_float|_optional_float]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:84-84`
- [[analysis.manifest_parsing__optional_float_tuple|_optional_float_tuple]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:293-293`
- [[analysis.manifest_parsing__parse_attitude_q|_parse_attitude_q]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:464-464`
- [[analysis.manifest_parsing__parse_vec3|_parse_vec3]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:455-455`
- [[analysis.manifest_parsing__require_float|_require_float]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:72-72`
- [[analysis.reporting__add_scaled_column_bang|_add_scaled_column!]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:59-59`
- [[analysis.reporting__append_display_error_columns_bang|_append_display_error_columns!]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:105-105`
- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:412-412`
- [[analysis.scenario_builders__make_tabulated_flight_density_model|_make_tabulated_flight_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:212-212`
- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:158-158`
- [[analysis.telemetry_loading__initial_time_et|_initial_time_et]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:30-30`
- [[analysis.telemetry_loading__load_time_aligned_telemetry|_load_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:382-382`
- [[analysis.telemetry_loading__sun_unit_vector_j2000|_sun_unit_vector_j2000]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:402-402`
- [[analysis.telemetry_loading__to_float_vector|_to_float_vector]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:9-9`
- [[assets.rpo_station_assets__stl_pointcloud|_stl_pointcloud]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:119-119`
- [[assets.rpo_station_assets_load_rpo_station_cad_triangles|load_rpo_station_cad_triangles]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:155-155`
- [[core.effector_sampling_statesample|StateSample]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:51-51`
- [[core.project_unit_quaternion|project_unit_quaternion]] · `callees` → `callers` · call · `src/core/numerics/quaternion_utils.jl:20-20`
- [[core.simulation_configuration_environmentmodel|EnvironmentModel]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:223-223`
- [[core.simulation_configuration_missionconfiguration|MissionConfiguration]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:166-166`
- [[dynamics.cloth_multibody__body_offset_velocity|_body_offset_velocity]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:183-183`
- [[dynamics.cloth_multibody_build_rectangular_compliant_grid|build_rectangular_compliant_grid]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:374-374`
- [[dynamics.cloth_multibody_rectangular_prism_inertia|rectangular_prism_inertia]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:189-189`
- [[dynamics.cloth_multibody_simulate_compliant_multibody|simulate_compliant_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:709-709`
- [[dynamics.cloth_multibody_step_compliant_multibody_implicit_midpoint|step_compliant_multibody_implicit_midpoint]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:655-655`
- [[dynamics.cloth_multibody_step_compliant_multibody_rk4|step_compliant_multibody_rk4]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:634-634`
- [[dynamics.cloth_robot_arm_dynamics__diag3|_diag3]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:98-98`
- [[dynamics.cloth_robot_arm_dynamics_cloth_reference_state|cloth_reference_state]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:32-32`
- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:466-466`
- [[dynamics.perturbations__harmonics_calcforcetorque_with_lpi|_harmonics_calcforcetorque_with_lpi]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1518-1518`
- [[dynamics.perturbations__harmonics_flat_batch_kernel_bang|_harmonics_flat_batch_kernel!]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:508-508`
- [[dynamics.perturbations__total_dipole_moment_body|_total_dipole_moment_body]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2062-2062`
- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:920-920`
- [[dynamics.perturbations_eddy_damping_torque|eddy_damping_torque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1917-1917`
- [[dynamics.perturbations_eddycurrentdampingmodel|EddyCurrentDampingModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1875-1875`
- [[dynamics.perturbations_lvlhcascadeattitudecontrolmodel|LVLHCascadeAttitudeControlModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:2163-2163`
- [[dynamics.perturbations_magnetictorquerodmodel|MagneticTorqueRodModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1991-1991`
- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:616-616`
- [[dynamics.point_mass_dynamics_mass_derivative|mass_derivative]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:19-19`
- [[dynamics.torque_models_body_angular_velocity|body_angular_velocity]] · `callees` → `callers` · call · `src/dynamics/rotational/torque_models.jl:6-6`
- [[dynamics.torque_models_body_torque|body_torque]] · `callees` → `callers` · call · `src/dynamics/rotational/torque_models.jl:2-2`
- [[dynx.rotational_attitude_kinematics_quaternion_derivative|quaternion_derivative]] · `callees` → `callers` · call · `src/dynamics/rotational/attitude_kinematics.jl:20-20`
- [[dynx.translational_point_mass_dynamics_acceleration_from_force|acceleration_from_force]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:5-5`
- [[dynx.translational_position_kinematics_position_derivative|position_derivative]] · `callees` → `callers` · call · `src/dynamics/translational/position_kinematics.jl:4-4`
- [[envana.ana_decay_diagnostics_secular_sma_slope|secular_sma_slope]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:64-64`
- [[environment.density_models__batch_elapsed_time|_batch_elapsed_time]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:840-840`
- [[environment.density_models__nrlmsise_ap_bins|_nrlmsise_ap_bins]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:534-534`
- [[environment.density_models__nrlmsise_ap_value|_nrlmsise_ap_value]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:320-320`
- [[environment.density_models__nrlmsise_eval_datetime|_nrlmsise_eval_datetime]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:510-510`
- [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:335-335`
- [[environment.density_models__planet_polyfit_valid_max_altitude_m|_planet_polyfit_valid_max_altitude_m]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:227-227`
- [[environment.density_models__planet_polyfit_valid_min_altitude_m|_planet_polyfit_valid_min_altitude_m]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:220-220`
- [[environment.density_models_exponentialatmospheremodel|ExponentialAtmosphereModel]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:56-56`
- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:938-938`
- [[environment.density_models_piecewiseexponentialatmospheremodel|PiecewiseExponentialAtmosphereModel]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:136-136`
- [[environment.density_models_polynomialfitatmospheremodel|PolynomialFitAtmosphereModel]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:210-210`
- [[environment.density_models_timetabulatedatmospheremodel|TimeTabulatedAtmosphereModel]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:773-773`
- [[environment.get_density|getDensity]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:773-773`
- [[environment.gravity_models__gravity_gradient_torque_body|_gravity_gradient_torque_body]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:87-87`
- [[environment.gravity_models__inverse_squared_gravity_accel|_inverse_squared_gravity_accel]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:40-40`
- [[environment.gravity_models__inverse_squared_j2_gravity_accel|_inverse_squared_j2_gravity_accel]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:46-46`
- [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:186-186`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:218-218`
- [[environment.simple_ephemerides__initial_time_datetime|_initial_time_datetime]] · `callees` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:49-49`
- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:138-138`
- [[ext.spaceagoragramsuiteext__gram_utc_string|_gram_utc_string]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:68-68`
- [[gnc.bridge_helpers__bridge_aerobraking_dry_mass|_bridge_aerobraking_dry_mass]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:110-110`
- [[gnc.bridge_helpers__bridge_aerobraking_entry_interface_m|_bridge_aerobraking_entry_interface_m]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:59-59`
- [[gnc.bridge_helpers__bridge_aerobraking_exit_interface_m|_bridge_aerobraking_exit_interface_m]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:68-68`
- [[gnc.bridge_helpers__bridge_aerobraking_max_heat_rate|_bridge_aerobraking_max_heat_rate]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:86-86`
- [[gnc.bridge_helpers__bridge_aerobraking_thrust_phi|_bridge_aerobraking_thrust_phi]] · `callees` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:122-122`
- [[gnc.clearance_rpo_path_clearance_stats|rpo_path_clearance_stats]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:18-18`
- [[gnc.config_robotarmsphereobstacle|RobotArmSphereObstacle]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/config.jl:8-8`
- [[gnc.constraint_tracking_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- [[gnc.control_commands_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:225-225`
- [[gnc.eom_predictor_f_ctrl_bang|f_ctrl!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:342-342`
- [[gnc.heat_load_control__edg_heat_load_scale_height|_edg_heat_load_scale_height]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:7-7`
- [[gnc.heat_load_control__edg_max_heat_load_for_links|_edg_max_heat_load_for_links]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:41-41`
- [[gnc.heat_load_control__edg_predict_mass|_edg_predict_mass]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:27-27`
- [[gnc.heat_load_control__edg_profile_heat_rates|_edg_profile_heat_rates]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:390-390`
- [[gnc.heat_load_control__edg_total_ref_area|_edg_total_ref_area]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:15-15`
- [[gnc.heat_load_control__edg_weighted_aero_coefficients|_edg_weighted_aero_coefficients]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:55-55`
- [[gnc.heat_rate_control__edg_maxwellian_heat_rate|_edg_maxwellian_heat_rate]] · `callees` → `callers` · call · `src/gnc/control/heat_rate_control.jl:110-110`
- [[gnc.hypr_utils_hypr_bezier_point|hypr_bezier_point]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:32-32`
- [[gnc.hypr_utils_hypr_iteration_weights|hypr_iteration_weights]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:96-96`
- [[gnc.hypr_utils_hypr_material_improvement|hypr_material_improvement]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:116-116`
- [[gnc.hypr_utils_hypr_protected_particle_mask|hypr_protected_particle_mask]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:133-133`
- [[gnc.hypr_utils_hypr_rrt_near_indices|hypr_rrt_near_indices]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:157-157`
- [[gnc.hypr_utils_hypr_rrt_steer|hypr_rrt_steer]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:169-169`
- [[gnc.lqmpc_rpo_discretize_zoh|rpo_discretize_zoh]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:30-30`
- [[gnc.lqmpc_rpo_hcw_continuous_mats|rpo_hcw_continuous_mats]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:3-3`
- [[gnc.lqmpc_rpo_ref_preview|rpo_ref_preview]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:130-130`
- [[gnc.path_costs_rpo_clearance_stats_from_samples|rpo_clearance_stats_from_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:43-43`
- [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:129-129`
- [[gnc.path_costs_rpo_obstacle_sigmoid_penalty|rpo_obstacle_sigmoid_penalty]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:76-76`
- [[gnc.path_retiming_rpo_interpolate_along_path|rpo_interpolate_along_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:42-42`
- [[gnc.path_retiming_rpo_remove_near_duplicate_samples|rpo_remove_near_duplicate_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:25-25`
- [[gnc.path_sampling_rpo_adaptive_sampling_min_ds_m|rpo_adaptive_sampling_min_ds_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:79-79`
- [[gnc.path_sampling_rpo_adaptive_sampling_step_m|rpo_adaptive_sampling_step_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:100-100`
- [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:130-130`
- [[gnc.path_sampling_rpo_inflated_obstacle_radius_m|rpo_inflated_obstacle_radius_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:69-69`
- [[gnc.path_sampling_rpo_sample_path_bezier|rpo_sample_path_bezier]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:25-25`
- [[gnc.path_sampling_rpo_sample_path_polyline|rpo_sample_path_polyline]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:61-61`
- [[gnc.planner_comparison__rpo_comparison_progress_bar|_rpo_comparison_progress_bar]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:69-69`
- [[gnc.planner_comparison_rpo_740_mpc_final_pso_config|rpo_740_mpc_final_pso_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:120-120`
- [[gnc.planner_comparison_rpo_comparison_cost_iteration_plot|rpo_comparison_cost_iteration_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1081-1081`
- [[gnc.planner_comparison_rpo_comparison_metric_value|rpo_comparison_metric_value]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:709-709`
- [[gnc.planner_comparison_rpo_group_metric_mean|rpo_group_metric_mean]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:682-682`
- [[gnc.planner_comparison_rpo_lqmpc_tracking_fuel_used_pct|rpo_lqmpc_tracking_fuel_used_pct]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:457-457`
- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:271-271`
- [[gnc.planner_core_robot_arm_hypr_path_cost_components|robot_arm_hypr_path_cost_components]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:19-19`
- [[gnc.propulsive_maneuvers__commanded_maneuver|_commanded_maneuver]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:103-103`
- [[gnc.propulsive_maneuvers__effective_direction_rad|_effective_direction_rad]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:169-169`
- [[gnc.propulsive_maneuvers__effective_thrust_isp|_effective_thrust_isp]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:180-180`
- [[gnc.propulsive_maneuvers__trace_maneuver_event_bang|_trace_maneuver_event!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:293-293`
- [[gnc.propulsive_maneuvers__validated_burn_plan|_validated_burn_plan]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:206-206`
- [[gnc.pso_adaptive_policy_rpo_estimate_geometry_complexity|rpo_estimate_geometry_complexity]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:7-7`
- [[gnc.pso_adaptive_policy_rpo_probe_geometry_metrics|rpo_probe_geometry_metrics]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:20-20`
- [[gnc.pso_parameters__rpo_pso_sync_sample_ds_with_safe_distance|_rpo_pso_sync_sample_ds_with_safe_distance]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:317-317`
- [[gnc.pso_parameters_rpo_hypr_refinement_sampling_density_m|rpo_hypr_refinement_sampling_density_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:332-332`
- [[gnc.pso_parameters_rpo_hypr_sampling_density_m|rpo_hypr_sampling_density_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:324-324`
- [[gnc.pso_path_planning_reset_swarm_bang|reset_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:323-323`
- [[gnc.pso_path_planning_rpo_pso_effective_safe_distance|rpo_pso_effective_safe_distance]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:106-106`
- [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:59-59`
- [[gnc.pso_refinement_rpo_refinement_segment_samples|rpo_refinement_segment_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:32-32`
- [[gnc.replanning__rpo_replanning_sphere|_rpo_replanning_sphere]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:71-71`
- [[gnc.replanning_rpo_active_replanning_spheres|rpo_active_replanning_spheres]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:111-111`
- [[gnc.replanning_rpo_plan_from_path|rpo_plan_from_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:257-257`
- [[gnc.replanning_rpo_replanning_sphere_center|rpo_replanning_sphere_center]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:106-106`
- [[gnc.replanning_rporeplanningsphere|RPOReplanningSphere]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:13-13`
- [[gnc.robot_arm_control_init_robot_arm_joint_mpc|init_robot_arm_joint_mpc]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:69-69`
- [[gnc.robot_arm_control_robot_arm_joint_mpc_reference_preview|robot_arm_joint_mpc_reference_preview]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:106-106`
- [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:146-146`
- [[gnc.rpo_guidance_hooks__rpo_record_replanning_event_bang|_rpo_record_replanning_event!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:71-71`
- [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:22-22`
- [[gnc.rrt_connect_rpo_rrt_collision_min_ds_m|rpo_rrt_collision_min_ds_m]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:75-75`
- [[gnc.rrt_connect_rpo_rrt_near_indices|rpo_rrt_near_indices]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:50-50`
- [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:100-100`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:505-505`
- [[gnc.rrt_warmstart__robot_arm_rrt_segment_samples|_robot_arm_rrt_segment_samples]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:35-35`
- [[gnc.station_geometry__rpo_build_station_kdtree|_rpo_build_station_kdtree]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl:38-38`
- [[gnc.struct_load_control__energy_depletion_struct_drag_area|_energy_depletion_struct_drag_area]] · `callees` → `callers` · call · `src/gnc/control/struct_load_control.jl:15-15`
- [[gnc.swarm_and_retiming__robot_arm_hypr_reaction_scale|_robot_arm_hypr_reaction_scale]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:137-137`
- [[gnc.swarm_and_retiming__robot_arm_hypr_reference_times_from_scales|_robot_arm_hypr_reference_times_from_scales]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:145-145`
- [[gnc.swarm_and_retiming__robot_arm_hypr_rigid_base_wrench_ratios|_robot_arm_hypr_rigid_base_wrench_ratios]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:183-183`
- [[gnc.target_energy_bracketing__edg_pos_vel_mass|_edg_pos_vel_mass]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:187-187`
- [[gnc.target_energy_bracketing_aerobrakingenergydepletionconfig|AerobrakingEnergyDepletionConfig]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:71-71`
- [[gnc.target_energy_bracketing_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:203-203`
- [[gnc.targeting_control__edg_certify_targeting_candidates|_edg_certify_targeting_candidates]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:891-891`
- [[gnc.targeting_control__edg_control_pos_vel_mass|_edg_control_pos_vel_mass]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:60-60`
- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:82-82`
- [[gnc.targeting_control__edg_in_drag_passage|_edg_in_drag_passage]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:92-92`
- [[gnc.targeting_control__edg_targeting_aero_acceleration|_edg_targeting_aero_acceleration]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:444-444`
- [[gnc.targeting_control__edg_targeting_constrained_alpha|_edg_targeting_constrained_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:382-382`
- [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:356-356`
- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:260-260`
- [[gnc.thruster_allocator_rpo_thruster_wrench_body|rpo_thruster_wrench_body]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/thruster_allocator.jl:24-24`
- [[gnc.thruster_guidance_functions__flight_apoapsis_ratio_scale|_flight_apoapsis_ratio_scale]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:108-108`
- [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:124-124`
- [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:212-212`
- [[gnc.trajectory_optimizers_rpo_chomp_numeric_gradient|rpo_chomp_numeric_gradient]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:212-212`
- [[gnc.trajectory_optimizers_rpo_chomp_obstacle_potential|rpo_chomp_obstacle_potential]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:159-159`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:379-379`
- [[gnc.trajectory_optimizers_rpo_trajectory_search_bounds|rpo_trajectory_search_bounds]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:112-112`
- [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:204-204`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:288-288`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:225-225`
- [[gncx.heat_rate_control__edg_heat_rate_alpha|_edg_heat_rate_alpha]] · `callees` → `callers` · call · `src/gnc/control/heat_rate_control.jl:143-143`
- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:516-516`
- [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl:19-19`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:147-147`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:186-186`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:325-325`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:245-245`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:288-288`
- [[gncz.clearance_robot_arm_clearance_stats_from_samples|robot_arm_clearance_stats_from_samples]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/clearance.jl:15-15`
- [[gncz.rpo_distance_queries_rpo_goal_standoff_point|rpo_goal_standoff_point]] · `callees` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/rpo_distance_queries.jl:6-6`
- [[gncz.rpo_plan_buffer_update_rpo_plan_buffer_bang|update_rpo_plan_buffer!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_plan_buffer.jl:23-23`
- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:227-227`
- [[grp.src_core_numerics|core/numerics/]] · `members_out` → `callers` · call · `src/core/numerics/quaternion_utils.jl:20-20`
- [[grp.src_core_state|core/state/]] · `members_out` → `callers` · call · `src/core/state/simulation_configuration.jl:223-223`
- [[grp.src_core_types|core/types/]] · `members_out` → `callers` · call · `src/core/types/effector_sampling.jl:51-51`
- [[grp.src_dynamics_coupled|dynamics/coupled/]] · `members_out` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1518-1518`
- [[grp.src_dynamics_multibody_cloth|dynamics/multibody_cloth/]] · `members_out` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:183-183`
- [[grp.src_dynamics_rotational|dynamics/rotational/]] · `members_out` → `callers` · call · `src/dynamics/rotational/torque_models.jl:6-6`
- [[grp.src_dynamics_translational|dynamics/translational/]] · `members_out` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:19-19`
- [[grp.src_environment_atmosphere|environment/atmosphere/]] · `members_out` → `callers` · call · `src/environment/atmosphere/density_models.jl:840-840`
- [[grp.src_environment_ephemerides|environment/ephemerides/]] · `members_out` → `callers` · call · `src/environment/ephemerides/simple_ephemerides.jl:49-49`
- [[grp.src_environment_gravity|environment/gravity/]] · `members_out` → `callers` · call · `src/environment/gravity/gravity_models.jl:87-87`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:205-205`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:238-238`
- [[grp.src_gnc_hypr|gnc/hypr/]] · `members_out` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:32-32`
- [[grp.src_gnc_internal|gnc/internal/]] · `members_out` → `callers` · call · `src/gnc/internal/bridge_helpers.jl:110-110`
- [[grp.src_gnc_navigation|gnc/navigation/]] · `members_out` → `callers` · call · `src/gnc/navigation/rpo_nav/distances/clearance.jl:18-18`
- [[grp.src_gnc_robotics|gnc/robotics/]] · `members_out` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/config.jl:8-8`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:151-151`
- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:75-75`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:178-178`
- [[io.io_serialization__load_checkpoint|_load_checkpoint]] · `callees` → `callers` · call · `src/io/serialization/io_serialization.jl:75-75`
- [[mission.maneuver_plans__phi_to_signed_maneuver_delta_v|_phi_to_signed_maneuver_delta_v]] · `callees` → `callers` · call · `src/mission/operations/maneuver_plans.jl:2-2`
- [[parallel.outer_route_selection__candidate_confidence_width|_candidate_confidence_width]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:472-472`
- [[parallel.outer_route_state__route_payload_stats|_route_payload_stats]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:119-119`
- [[parallel.outer_route_state__route_stats_payload|_route_stats_payload]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:96-96`
- [[parallel.persistent_hints__hint_mean_and_width|_hint_mean_and_width]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:235-235`
- [[parallel.persistent_hints__hint_payload_stats|_hint_payload_stats]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:37-37`
- [[parallel.persistent_hints__hint_record_observation_bang|_hint_record_observation!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:352-352`
- [[parallel.persistent_hints__hint_stats_payload|_hint_stats_payload]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:14-14`
- [[parcore.no_gram_presets_make_no_gram_environment|make_no_gram_environment]] · `callees` → `callers` · call · `src/core/state/no_gram_presets.jl:90-90`
- [[simulation.campaign_route_features|campaign_route_features]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:75-75`
- [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:178-178`
- [[simulation.dynamics_rhs__rhs_effector_static_cost_ns|_rhs_effector_static_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:352-352`
- [[simulation.dynamics_rhs__robot_arm_effector_matches|_robot_arm_effector_matches]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1629-1629`
- [[simulation.dynamics_rhs__update_rhs_flat_packet_cost_model_bang|_update_rhs_flat_packet_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:633-633`
- [[simulation.dynamics_rhs__update_rhs_flat_packet_overhead_model_bang|_update_rhs_flat_packet_overhead_model!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:600-600`
- [[simulation.effector_sampling__extract_sample_mass_kg|_extract_sample_mass_kg]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:16-16`
- [[simulation.event_callbacks_affect_upcrossing_bang|affect_upcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:151-151`
- [[simulation.execution__append_backbone_saved_segment_bang|_append_backbone_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:77-77`
- [[simulation.execution__append_checkpoint_saved_segment_bang|_append_checkpoint_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:95-95`
- [[simulation.model_selection__gram_batch_elapsed_time|_gram_batch_elapsed_time]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:83-83`
- [[simulation.monte_carlo_result|MonteCarloResult]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:62-62`
- [[simulation.refresh__gram_track_cache_fill_from_trajectory_bang|_gram_track_cache_fill_from_trajectory!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:23-23`
- [[simulation.rhs_calibration__rhs_calib_load_bang|_rhs_calib_load!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:121-121`
- [[simulation.rhs_calibration__rhs_calib_save_bang|_rhs_calib_save!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:139-139`
- [[simulation.rhs_calibration__run_rhs_sweep_bang|_run_rhs_sweep!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:282-282`
- [[simulation.runtime__density_segment_end_t|_density_segment_end_t]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:49-49`
- [[simulation.runtime__extract_mass_kg|_extract_mass_kg]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:7-7`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:293-293`
- [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`
- [[simulation.save_fields__save_altitude|_save_altitude]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:87-87`
- [[simulation.save_fields__save_heat_rate|_save_heat_rate]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:132-132`
- [[simulation.save_fields__save_latitude_deg|_save_latitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:101-101`
- [[simulation.save_fields__save_longitude_deg|_save_longitude_deg]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:115-115`
- [[simulation.setup__cache_from_nbody_ephemeris_payload|_cache_from_nbody_ephemeris_payload]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1641-1641`
- [[simulation.setup__effector_observed_cost_ns_per_item|_effector_observed_cost_ns_per_item]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:520-520`
- [[simulation.setup__ephemerides_time_seconds_flexible|_ephemerides_time_seconds_flexible]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1507-1507`
- [[simulation.setup__nbody_ephemeris_cache_payload|_nbody_ephemeris_cache_payload]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1604-1604`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1692-1692`
- [[simulation.setup__rhs_effector_observed_cost_ns|_rhs_effector_observed_cost_ns]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:593-593`
- [[simulation.setup__update_effector_cost_model_bang|_update_effector_cost_model!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:543-543`
- [[simulation.solver_policy__solve_with_explicit_solver|_solve_with_explicit_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:341-341`
- [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:549-549`
- [[simulation.solver_policy__solve_with_multirate_solver|_solve_with_multirate_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:411-411`
- [[simulation.state_access__state_mass_kg|_state_mass_kg]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:65-65`
- [[simulation.targeting__gram_entry_mass_kg|_gram_entry_mass_kg]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:211-211`
- [[simulation.targeting__gram_entry_reference_area_m2|_gram_entry_reference_area_m2]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:190-190`
- [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:81-81`
- [[simulation.targeting__gram_orbit_period_target|_gram_orbit_period_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:152-152`
- [[simulation.targeting__gram_periapsis_target|_gram_periapsis_target]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:120-120`
- [[simulation.thermal_callbacks_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:76-76`
- [[simulation.vacuum_predicted_gram__vacuum_j2_accel|_vacuum_j2_accel]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:65-65`
- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:204-204`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:76-76`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:385-385`
- [[vehicle.model__initial_condition_radius_for_oblate_altitude|_initial_condition_radius_for_oblate_altitude]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:198-198`
- [[vehicle.model_initialcondition|InitialCondition]] · `callees` → `callers` · call · `src/vehicle/spacecraft/model.jl:53-53`
- [[vehicle.robotics_closest_surface_target|closest_surface_target]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:315-315`
- [[vehicle.robotics_default_cloth_arm_model|default_cloth_arm_model]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:131-131`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The conversion relies on the CSV columns `telemetry_axis`, `sim_interp_value_km` and `error_km` existing and being numeric; a missing column raises a `DataFrame` field error rather than a descriptive message. Using a `Float64` time value as a `Dict` key means sample matching across runs requires bit-identical times; any floating drift in `telemetry_axis` between runs silently drops that sample from `axes_common`, and fewer than 12 common samples aborts the fit.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/ic_fit.jl` line 24.
