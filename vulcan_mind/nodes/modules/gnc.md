---
id: module.gnc
label: CommandTypes
kind: module
source:
  file: src/gnc/command_types.jl
  symbol: CommandTypes
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: Exported guidance, navigation and control surface — command records
    `PropulsiveManeuverCommand`, `PropulsiveBurnPlan` and `AerobrakingControlCommand`,
    the RPO controller entry points `init_rpo_lqmpc` and `rpo_lqmpc_control`, and
    the arm planners `plan_robot_arm_motion_hypr` and `cloth_ik` consumers.
tags:
- module
charts:
- master
origin: agent
---

# CommandTypes

## Purpose
`CommandTypes` is the command-record root of the SpaceAGORA GNC subsystem under `src/gnc/`. It declares the plain data structures that guidance produces and control consumes, so that guidance solvers, control effectors and the ODE right-hand side never share mutable state. The wider module tree it anchors spans aerobraking guidance (`src/gnc/guidance/aerobraking/`), rendezvous and proximity operations model-predictive control (`src/gnc/control/rpo_mpc/`), robot-arm planning (`src/gnc/robotics/`) and relative navigation (`src/gnc/navigation/rpo_nav`).

## Theory & Math
The dominant control law in this subsystem is the finite-horizon linear-quadratic model predictive controller in `src/gnc/control/rpo_mpc/lqmpc.jl`. Relative translation is modelled with the Clohessy-Wiltshire equations, written in the target's radial-transverse-normal frame:

$$
\ddot{x} = 3n^2 x + 2n\dot{y} + u_x,\qquad
\ddot{y} = -2n\dot{x} + u_y,\qquad
\ddot{z} = -n^2 z + u_z
$$

where $x$ is the radial offset in metres, $y$ the along-track offset in metres, $z$ the cross-track offset in metres, $n$ the target mean motion in rad/s, and $u = (u_x,u_y,u_z)$ the commanded chaser acceleration in m/s². Stacking $\mathbf{x} = [x,y,z,\dot{x},\dot{y},\dot{z}]^\top$ gives $\dot{\mathbf{x}} = A\mathbf{x} + Bu$, exactly the `A` and `B` built by `rpo_hcw_continuous_mats`.

The controller minimises

$$
J = \sum_{k=1}^{N-1} (\mathbf{x}_k - \mathbf{r}_k)^\top Q (\mathbf{x}_k - \mathbf{r}_k)
  + (\mathbf{x}_N - \mathbf{r}_N)^\top Q_f (\mathbf{x}_N - \mathbf{r}_N)
  + \sum_{k=0}^{N-1} u_k^\top R\, u_k
$$

subject to $\mathbf{x}_{k+1} = A_d\mathbf{x}_k + B_d u_k$ and $u_{\min} \le u_k \le u_{\max}$, where $N$ is `horizon`, $\mathbf{r}_k$ the reference state produced by `rpo_ref_preview`, $Q$ and $Q_f$ the stage and terminal state weights, and $R$ the control weight.

The robot-arm branch solves the inverse geometry problem by damped least squares, $\Delta q = (J^\top J + \lambda^2 I)^{-1} J^\top e$, with $J$ the end-effector position Jacobian, $e$ the Cartesian position error in metres and $\lambda$ the damping factor.

## Model & Assumptions
- Clohessy-Wiltshire linearisation assumes a circular target orbit and a chaser separation small relative to the target orbit radius; eccentricity and separations of many kilometres break the $3n^2x$ and $-n^2z$ gravity-gradient terms.
- Control acceleration is treated as a zero-order hold over `dt`, matched by `rpo_discretize_zoh`; commands held longer than one control tick accumulate prediction error.
- `PropulsiveManeuverCommand` and `PropulsiveBurnPlan` are immutable `Base.@kwdef` structs with `valid = false` defaults, so a guidance model that fails to converge returns an inert command rather than an exception.
- `AerobrakingControlCommand` carries a single scalar `alpha_command` in radians, which assumes attitude control resolves angle of attack without a separate bank-angle channel.

## Design & Implementation
`src/gnc/command_types.jl` (lines 1–29) defines three keyword-constructed immutable structs and exports them on line 3. `PropulsiveBurnPlan` extends `PropulsiveManeuverCommand` with burn-window timing (`start_burn_s`, `stop_burn_s`), engine parameters (`thrust_n`, `isp_s`) and the derived `commanded_impulse_n_s` and `propellant_required_kg`, so the propagator can integrate mass depletion without re-deriving the plan.

`src/gnc/control/control_hooks.jl` is the assembly point: it imports `CommandTypes` on line 10, then `include`s `rpo_mpc/lqmpc.jl`, `rpo_mpc/rpo_mpc_control_model.jl`, `robot_arm_control.jl`, `momentum_manager.jl` and the three aerobraking files. `calcControlEffect!` in `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl` (lines 7–36) is the per-step driver: it forms the relative state with `inertial_to_rtn_relative_state`, previews the reference with `rpo_ref_preview`, calls `rpo_lqmpc_control`, rotates the resulting acceleration into the body frame, and allocates it across six thrusters with `rpo_allocate_six_axis_thrusters`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | Exported guidance, navigation and control surface — command records `PropulsiveManeuverCommand`, `PropulsiveBurnPlan` and `AerobrakingControlCommand`, the RPO controller entry points `init_rpo_lqmpc` and `rpo_lqmpc_control`, and the arm planners `plan_robot_arm_motion_hypr` and `cloth_ik` consumers. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[gnc.bridge_helpers__bridge_optional_field|_bridge_optional_field]] · `module_api` · call · `src/gnc/internal/bridge_helpers.jl`
- `api` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `module_api` · call · `src/gnc/internal/bridge_helpers.jl`
- `api` → [[gnc.command_types_aerobrakingcontrolcommand|AerobrakingControlCommand]] · `module_api` · call · `src/gnc/command_types.jl`
- `api` → [[gnc.command_types_commandtypes|CommandTypes]] · `module_api` · call · `src/gnc/command_types.jl`
- `api` → [[gnc.config_robotarmsphereobstacle|RobotArmSphereObstacle]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/config.jl`
- `api` → [[gnc.control_commands_asim_ctrl_rf|asim_ctrl_rf]] · `module_api` · call · `src/gnc/control/aerobraking/control_commands.jl`
- `api` → [[gnc.heat_load_control__edg_closed_form_heat_load_trajectory|_edg_closed_form_heat_load_trajectory]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_eccentric_anomaly_from_true|_edg_eccentric_anomaly_from_true]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_first_low_alpha_interval_indices|_edg_first_low_alpha_interval_indices]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_heat_load_alpha_profile|_edg_heat_load_alpha_profile]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_heat_load_lambdas|_edg_heat_load_lambdas]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_heat_load_low_drag_active|_edg_heat_load_low_drag_active]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_heat_load_scale_height|_edg_heat_load_scale_height]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_heat_load_track_env|_edg_heat_load_track_env]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_integrate_series|_edg_integrate_series]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_integrated_heat_load_trajectory|_edg_integrated_heat_load_trajectory]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_max_heat_load_for_links|_edg_max_heat_load_for_links]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_mean_anomaly_from_true|_edg_mean_anomaly_from_true]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_prediction_time_grid|_edg_prediction_time_grid]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_profile_heat_rates|_edg_profile_heat_rates]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_sample_prediction_atmosphere|_edg_sample_prediction_atmosphere]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_total_ref_area|_edg_total_ref_area]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control__edg_weighted_aero_coefficients|_edg_weighted_aero_coefficients]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_load_control_acceleration|acceleration]] · `module_api` · call · `src/gnc/control/heat_load_control.jl`
- `api` → [[gnc.heat_rate_control__edg_maxwellian_heat_rate|_edg_maxwellian_heat_rate]] · `module_api` · call · `src/gnc/control/heat_rate_control.jl`
- `api` → [[gnc.heat_rate_control__energy_depletion_heat_rate_calc|_energy_depletion_heat_rate_calc]] · `module_api` · call · `src/gnc/control/heat_rate_control.jl`
- `api` → [[gnc.heat_rate_models_aoa|aoa]] · `module_api` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl`
- `api` → [[gnc.hypr_utils__hypr_rrt_costs|_hypr_rrt_costs]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils__hypr_rrt_parents|_hypr_rrt_parents]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_bezier_point|hypr_bezier_point]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_bezier_point_bang|hypr_bezier_point!]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_iteration_weights|hypr_iteration_weights]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_material_improvement|hypr_material_improvement]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_path_length|hypr_path_length]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_protected_particle_mask|hypr_protected_particle_mask]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_rrt_join_paths|hypr_rrt_join_paths]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_rrt_near_indices|hypr_rrt_near_indices]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_rrt_nearest_index|hypr_rrt_nearest_index]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_rrt_refresh_subtree_costs_bang|hypr_rrt_refresh_subtree_costs!]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_rrt_steer|hypr_rrt_steer]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hypr_rrt_tree_path|hypr_rrt_tree_path]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.hypr_utils_hyprutils|HYPRUtils]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gnc.interfaces_abstractaerobrakingstrategy|AbstractAerobrakingStrategy]] · `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`
- `api` → [[gnc.interfaces_aerobrakingguidanceinput|AerobrakingGuidanceInput]] · `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`
- `api` → [[gnc.interfaces_tedgstrategy|TEdgStrategy]] · `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`
- `api` → [[gnc.lqmpc_init_rpo_lqmpc|init_rpo_lqmpc]] · `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`
- `api` → [[gnc.lqmpc_rpo_block_diag|rpo_block_diag]] · `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`
- `api` → [[gnc.lqmpc_rpo_discretize_zoh|rpo_discretize_zoh]] · `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`
- `api` → [[gnc.lqmpc_rpo_hcw_continuous_mats|rpo_hcw_continuous_mats]] · `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`
- `api` → [[gnc.lqmpc_rpo_prediction_mats|rpo_prediction_mats]] · `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`
- `api` → [[gnc.lqmpc_rpolqmpccontroller|RpoLQMPCController]] · `module_api` · call · `src/gnc/control/rpo_mpc/lqmpc.jl`
- `api` → [[gnc.mesh_distance__nearest_station_distance_sq_kdtree|_nearest_station_distance_sq_kdtree]] · `module_api` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl`
- `api` → [[gnc.mesh_distance_nearest_station_distance_sq|nearest_station_distance_sq]] · `module_api` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl`
- `api` → [[gnc.momentum_manager_calccontrolforcetorque|calcControlForceTorque]] · `module_api` · call · `src/gnc/control/momentum_manager.jl`
- `api` → [[gnc.navigation_hooks_navigationhooks|NavigationHooks]] · `module_api` · call · `src/gnc/navigation/navigation_hooks.jl`
- `api` → [[gnc.path_costs_rpo_clearance_stats_from_samples|rpo_clearance_stats_from_samples]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`
- `api` → [[gnc.path_costs_rpo_fuel_proxy_from_samples|rpo_fuel_proxy_from_samples]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`
- `api` → [[gnc.path_costs_rpo_obstacle_sigmoid_penalty|rpo_obstacle_sigmoid_penalty]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`
- `api` → [[gnc.path_costs_rpo_path_cost_normalization_refs|rpo_path_cost_normalization_refs]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl`
- `api` → [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`
- `api` → [[gnc.path_sampling_rpo_bezier_point|rpo_bezier_point]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`
- `api` → [[gnc.path_sampling_rpo_inflated_obstacle_radius_m|rpo_inflated_obstacle_radius_m]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`
- `api` → [[gnc.path_sampling_rpo_sample_path_bezier|rpo_sample_path_bezier]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`
- `api` → [[gnc.path_sampling_rpo_sample_path_polyline_adaptive|rpo_sample_path_polyline_adaptive]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl`
- `api` → [[gnc.planner_comparison__rpo_comparison_progress_bar|_rpo_comparison_progress_bar]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison__rpo_comparison_rrt_connect_seed_path|_rpo_comparison_rrt_connect_seed_path]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison__rpo_comparison_rrt_connect_settings|_rpo_comparison_rrt_connect_settings]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison__rpo_comparison_rrt_star_settings|_rpo_comparison_rrt_star_settings]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_740_mpc_final_pso_config|rpo_740_mpc_final_pso_config]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_artifact_slug|rpo_comparison_artifact_slug]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_axis_name|rpo_comparison_axis_name]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_cost_iteration_plot|rpo_comparison_cost_iteration_plot]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_cost_iteration_xvalues|rpo_comparison_cost_iteration_xvalues]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_failed_paths_plot|rpo_comparison_failed_paths_plot]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_metric_includes_result|rpo_comparison_metric_includes_result]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_metric_specs|rpo_comparison_metric_specs]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_metric_value|rpo_comparison_metric_value]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_path_family_plot|rpo_comparison_path_family_plot]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_single_path_plot|rpo_comparison_single_path_plot]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_station_mesh_trace|rpo_comparison_station_mesh_trace]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_station_trace|rpo_comparison_station_trace]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_comparison_trace_axis_name|rpo_comparison_trace_axis_name]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_flatten_planner_results|rpo_flatten_planner_results]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_group_metric_mean|rpo_group_metric_mean]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_lqmpc_reference_preview|rpo_lqmpc_reference_preview]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_lqmpc_tracking_fuel_used_pct|rpo_lqmpc_tracking_fuel_used_pct]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_write_failed_path_outputs|rpo_write_failed_path_outputs]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_write_namedtuple_csv|rpo_write_namedtuple_csv]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpo_write_planner_comparison_outputs|rpo_write_planner_comparison_outputs]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpolqmpctrackingsettings|RPOLQMPCTrackingSettings]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_comparison_rpoplannercomparisoncase|RPOPlannerComparisonCase]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`
- `api` → [[gnc.planner_core__robot_arm_hypr_refinement_better|_robot_arm_hypr_refinement_better]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl`
- `api` → [[gnc.propulsive_maneuver_command|PropulsiveManeuverCommand]] · `module_api` · call · `src/gnc/command_types.jl`
- `api` → [[gnc.propulsive_maneuvers__available_propellant_kg|_available_propellant_kg]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__burn_plan_buffer|_burn_plan_buffer]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__control_effector_log_enabled|_control_effector_log_enabled]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__control_effector_strict_exceptions|_control_effector_strict_exceptions]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__effective_burn_window|_effective_burn_window]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__effective_direction_rad|_effective_direction_rad]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__effective_thrust_isp|_effective_thrust_isp]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__guidance_maneuver_command|_guidance_maneuver_command]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__maneuver_trace_enabled|_maneuver_trace_enabled]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__maneuver_trace_path|_maneuver_trace_path]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__safe_orbit_counter|_safe_orbit_counter]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers__trace_bool_enabled|_trace_bool_enabled]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers_calccontrolmassflowrate|calcControlMassFlowRate]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.propulsive_maneuvers_calcreactionwheeltorque|calcReactionWheelTorque]] · `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`
- `api` → [[gnc.pso_adaptive_policy_rpo_probe_geometry_metrics|rpo_probe_geometry_metrics]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl`
- `api` → [[gnc.pso_parameters__rpo_pso_config_tuple|_rpo_pso_config_tuple]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters__rpo_pso_normalize_kwargs|_rpo_pso_normalize_kwargs]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters__rpo_pso_sync_sample_ds_with_safe_distance|_rpo_pso_sync_sample_ds_with_safe_distance]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpo_hypr_refinement_sampling_density_m|rpo_hypr_refinement_sampling_density_m]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpo_hypr_sampling_density_m|rpo_hypr_sampling_density_m]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpoadaptivesamplingsettings|RPOAdaptiveSamplingSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoadaptivesettings|RPOPSOAdaptiveSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsocullsettings|RPOPSOCullSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoearlystoppingsettings|RPOPSOEarlyStoppingSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoobjectivesettings|RPOPSOObjectiveSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoprobesettings|RPOPSOProbeSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoreexploresettings|RPOPSOReexploreSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsorefinementsettings|RPOPSORefinementSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoretimingsettings|RPOPSORetimingSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsorrtconnectwarmstartsettings|RPOPSORRTConnectWarmstartSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoschedulesettings|RPOPSOScheduleSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsostagnationsettings|RPOPSOStagnationSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_parameters_rpopsoswarmsettings|RPOPSOSwarmSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`
- `api` → [[gnc.pso_path_planning_rpo_pso_project_to_segment|rpo_pso_project_to_segment]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl`
- `api` → [[gnc.pso_path_planning_rpo_pso_tapered_noise_scale|rpo_pso_tapered_noise_scale]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_bernstein|rpo_refinement_bernstein]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_better|rpo_refinement_better]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_clamp_path|rpo_refinement_clamp_path]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_project_to_segment|rpo_refinement_project_to_segment]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_sample_params|rpo_refinement_sample_params]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_segment_samples|rpo_refinement_segment_samples]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_refinement_shortcut_samples|rpo_refinement_shortcut_samples]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.pso_refinement_rpo_try_accept_refinement|rpo_try_accept_refinement]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`
- `api` → [[gnc.replanning__rpo_replanning_property|_rpo_replanning_property]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning__rpo_replanning_sphere|_rpo_replanning_sphere]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning_rpo_plan_from_path|rpo_plan_from_path]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning_rpo_replanning_signature|rpo_replanning_signature]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning_rpo_replanning_sphere_center|rpo_replanning_sphere_center]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning_rpo_retime_existing_plan|rpo_retime_existing_plan]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning_rpo_sphere_surface_points|rpo_sphere_surface_points]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning_rporeplanningconfig|RPOReplanningConfig]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.replanning_rporeplanningsphere|RPOReplanningSphere]] · `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`
- `api` → [[gnc.robot_arm_control__robot_arm_control_axis_angle_about|_robot_arm_control_axis_angle_about]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control__robot_arm_control_quat_conj|_robot_arm_control_quat_conj]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control__robot_arm_control_reference_state|_robot_arm_control_reference_state]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control__robot_arm_control_spacecraft_state|_robot_arm_control_spacecraft_state]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control__robot_arm_default_joint_inertia|_robot_arm_default_joint_inertia]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_calccontrolforcetorque|calcControlForceTorque]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_calccontrolmassflowrate|calcControlMassFlowRate]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_init_robot_arm_joint_mpc|init_robot_arm_joint_mpc]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_robot_arm_joint_mpc_reference_preview|robot_arm_joint_mpc_reference_preview]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_robot_arm_measured_joint_state|robot_arm_measured_joint_state]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_robotarmcontroleffector|RobotArmControlEffector]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_robotarmheldactuation|RobotArmHeldActuation]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_control_robotarmjointmpccontroller|RobotArmJointMPCController]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gnc.robot_arm_planning_robot_arm_plan_sample|robot_arm_plan_sample]] · `module_api` · call · `src/gnc/robotics/robot_arm_planning.jl`
- `api` → [[gnc.robot_arm_planning_robotarmplanning|RobotArmPlanning]] · `module_api` · call · `src/gnc/robotics/robot_arm_planning.jl`
- `api` → [[gnc.rpo_guidance_hooks__rpo_record_replanning_event_bang|_rpo_record_replanning_event!]] · `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`
- `api` → [[gnc.rpo_guidance_hooks__rpo_replanning_config|_rpo_replanning_config]] · `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`
- `api` → [[gnc.rpo_guidance_hooks__rpo_state_pos_vel|_rpo_state_pos_vel]] · `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`
- `api` → [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`
- `api` → [[gnc.rpo_mpc_control_model_calccontrolforcetorque|calcControlForceTorque]] · `module_api` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl`
- `api` → [[gnc.rpo_mpc_control_model_calccontrolmassflowrate|calcControlMassFlowRate]] · `module_api` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl`
- `api` → [[gnc.rpo_mpc_control_model_calcreactionwheeltorque|calcReactionWheelTorque]] · `module_api` · call · `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl`
- `api` → [[gnc.rpo_plan_buffer_rpoplan|RPOPlan]] · `module_api` · call · `src/gnc/guidance/rpo/rpo_plan_buffer.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_collision_min_ds_m|rpo_rrt_collision_min_ds_m]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_near_indices|rpo_rrt_near_indices]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_nearest_index|rpo_rrt_nearest_index]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_refresh_subtree_costs_bang|rpo_rrt_refresh_subtree_costs!]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_star_add_node_bang|rpo_rrt_star_add_node!]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_steer|rpo_rrt_steer]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rpo_rrt_tree_path|rpo_rrt_tree_path]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_connect_rporrtstarsettings|RPORRTStarSettings]] · `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`
- `api` → [[gnc.rrt_warmstart__robot_arm_rrt_nearest_index|_robot_arm_rrt_nearest_index]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl`
- `api` → [[gnc.rrt_warmstart__robot_arm_rrt_segment_samples|_robot_arm_rrt_segment_samples]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl`
- `api` → [[gnc.rrt_warmstart__robot_arm_rrt_steer|_robot_arm_rrt_steer]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl`
- `api` → [[gnc.station_geometry__rpo_build_station_kdtree|_rpo_build_station_kdtree]] · `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl`
- `api` → [[gnc.station_geometry_rpostationkdnode|RPOStationKDNode]] · `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl`
- `api` → [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_base_wrench_ratios|_robot_arm_hypr_cloth_base_wrench_ratios]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`
- `api` → [[gnc.swarm_and_retiming__robot_arm_hypr_cloth_state_for_reaction|_robot_arm_hypr_cloth_state_for_reaction]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`
- `api` → [[gnc.swarm_and_retiming__robot_arm_hypr_link_com_history|_robot_arm_hypr_link_com_history]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`
- `api` → [[gnc.swarm_and_retiming__robot_arm_hypr_rigid_base_wrench_ratios|_robot_arm_hypr_rigid_base_wrench_ratios]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`
- `api` → [[gnc.swarm_and_retiming__robot_arm_path_length|_robot_arm_path_length]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`
- `api` → [[gnc.swarm_and_retiming__robot_arm_path_smoothness|_robot_arm_path_smoothness]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl`
- `api` → [[gnc.target_energy_bracketing__edg_interpolate_bracket_value|_edg_interpolate_bracket_value]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.target_energy_bracketing__edg_state_index_ok|_edg_state_index_ok]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.target_energy_bracketing__edg_symbol_tuple|_edg_symbol_tuple]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.target_energy_bracketing__edg_validate_symbol_set|_edg_validate_symbol_set]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.target_energy_bracketing_aerobrakingenergydepletionconfig|AerobrakingEnergyDepletionConfig]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.target_energy_bracketing_aerobrakingenergydepletionguidancemodel|AerobrakingEnergyDepletionGuidanceModel]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.target_energy_bracketing_aerobrakingenergydepletionstate|AerobrakingEnergyDepletionState]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.target_energy_bracketing_calcguidanceeffect_bang|calcGuidanceEffect!]] · `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`
- `api` → [[gnc.targeting_control__apply_solar_panel_aoa_bang|_apply_solar_panel_aoa!]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_base_alpha|_edg_base_alpha]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_certify_targeting_candidates|_edg_certify_targeting_candidates]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_command_alpha_bang|_edg_command_alpha!]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_control_pos_vel_mass|_edg_control_pos_vel_mass]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_control_sat_state|_edg_control_sat_state]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_control_state_index_ok|_edg_control_state_index_ok]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_disable_uncertified_targeting_bang|_edg_disable_uncertified_targeting!]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_ephemeris_time|_edg_ephemeris_time]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_integrated_max_energy_depletion_trajectory|_edg_integrated_max_energy_depletion_trajectory]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_integrated_targeting_trajectory|_edg_integrated_targeting_trajectory]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_orbit_metrics_from_rv|_edg_orbit_metrics_from_rv]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_panel_link_tuple|_edg_panel_link_tuple]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_planet_frame_lpi|_edg_planet_frame_lpi]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_recompute_switches_bang|_edg_recompute_switches!]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_solve_targeting_switch|_edg_solve_targeting_switch]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_target_energy_from_apoapsis|_edg_target_energy_from_apoapsis]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_targeting_aero_acceleration|_edg_targeting_aero_acceleration]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_targeting_bracket_outcomes|_edg_targeting_bracket_outcomes]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_targeting_constrained_alpha|_edg_targeting_constrained_alpha]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_targeting_prediction_time_grid|_edg_targeting_prediction_time_grid]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control__edg_targeting_switch_outcomes|_edg_targeting_switch_outcomes]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_acceleration|acceleration]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_apoapsis_residual|apoapsis_residual]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_calccontrolforcetorque|calcControlForceTorque]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_energy_residual|energy_residual]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_max_energy_alpha|max_energy_alpha]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_solarpanelangleofattackcontrolmodel|SolarPanelAngleOfAttackControlModel]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_solve_apoapsis_switch|solve_apoapsis_switch]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_control_solve_energy_switch|solve_energy_switch]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gnc.targeting_solver_control_solarpanels_targeting_closed_form|control_solarpanels_targeting_closed_form]] · `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`
- `api` → [[gnc.targeting_solver_control_solarpanels_targeting_heatload|control_solarpanels_targeting_heatload]] · `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`
- `api` → [[gnc.targeting_solver_control_solarpanels_targeting_num_int|control_solarpanels_targeting_num_int]] · `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`
- `api` → [[gnc.targeting_solver_func_e|func_e]] · `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`
- `api` → [[gnc.targeting_solver_func_targeting_heatload|func_targeting_heatload]] · `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`
- `api` → [[gnc.targeting_solver_func_targeting_num_int|func_targeting_num_int]] · `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`
- `api` → [[gnc.thruster_guidance_functions__ensure_apo_target_state_bang|_ensure_apo_target_state!]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gnc.thruster_guidance_functions__flight_apoapsis_ratio_scale|_flight_apoapsis_ratio_scale]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gnc.thruster_guidance_functions__oblate_altitude_from_radius|_oblate_altitude_from_radius]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gnc.thruster_guidance_functions__oblate_surface_radius|_oblate_surface_radius]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gnc.thruster_guidance_functions__radius_for_oblate_altitude|_radius_for_oblate_altitude]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gnc.thruster_guidance_functions__wrap_2pi_guidance|_wrap_2pi_guidance]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gnc.thruster_guidance_models_apoapsistargetperiapsisraiseguidancemodel|ApoapsisTargetPeriapsisRaiseGuidanceModel]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl`
- `api` → [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`
- `api` → [[gnc.tracking_executor__control_solarpanels_openloop_impl|_control_solarpanels_openloop_impl]] · `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`
- `api` → [[gnc.tracking_executor__control_strict_exceptions|_control_strict_exceptions]] · `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`
- `api` → [[gnc.tracking_executor_control_solarpanels_heatload|control_solarpanels_heatload]] · `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`
- `api` → [[gnc.tracking_executor_control_solarpanels_openloop|control_solarpanels_openloop]] · `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`
- `api` → [[gnc.tracking_executor_no_control|no_control]] · `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`
- `api` → [[gnc.trajectory_optimizers_rpo_chomp_obstacle_potential|rpo_chomp_obstacle_potential]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`
- `api` → [[gnc.trajectory_optimizers_rpo_soft_obstacle_cost_from_samples|rpo_soft_obstacle_cost_from_samples]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`
- `api` → [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`
- `api` → [[gnc.trajectory_optimizers_rpo_stomp_waypoint_state_cost|rpo_stomp_waypoint_state_cost]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`
- `api` → [[gnc.trajectory_optimizers_rpo_trajectory_internal_points|rpo_trajectory_internal_points]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`
- `api` → [[gnc.trajectory_optimizers_rpo_trajectory_smoothness_cost|rpo_trajectory_smoothness_cost]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`
- `api` → [[gnc.trajectory_optimizers_rpostompsettings|RPOSTOMPSettings]] · `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`
- `api` → [[gnc.trajectory_predictor_closed_form_targeting|closed_form_targeting]] · `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl`
- `api` → [[gncx.control_hooks_controlhooks|ControlHooks]] · `module_api` · call · `src/gnc/control/control_hooks.jl`
- `api` → [[gncx.momentum_manager_magneticmomentummanagermodel|MagneticMomentumManagerModel]] · `module_api` · call · `src/gnc/control/momentum_manager.jl`
- `api` → [[gncx.robot_arm_control_robot_arm_joint_mpc_control|robot_arm_joint_mpc_control]] · `module_api` · call · `src/gnc/control/robot_arm_control.jl`
- `api` → [[gncx.targeting_control_aerobrakingenergydepletioncontrolmodel|AerobrakingEnergyDepletionControlModel]] · `module_api` · call · `src/gnc/control/targeting_control.jl`
- `api` → [[gncy.guidance_models_guidancemodels|GuidanceModels]] · `module_api` · call · `src/gnc/guidance/guidance_models.jl`
- `api` → [[gncy.interfaces_eedgstrategy|EEdgStrategy]] · `module_api` · call · `src/gnc/guidance/aerobraking/interfaces.jl`
- `api` → [[gncz.cubesat_geometry_rpocubesatgeometry|RPOCubeSatGeometry]] · `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/cubesat_geometry.jl`
- `api` → [[gncz.hypr_utils_hypr_sample_count_path|hypr_sample_count_path]] · `module_api` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[gncz.robot_arm_hypr_robot_arm_hypr|robot_arm_hypr]] · `module_api` · call · `src/gnc/robotics/robot_arm_hypr.jl`
- `api` → [[gncz.rpo_reference_geometry_rporeferencegeometry|RPOReferenceGeometry]] · `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/rpo_reference_geometry.jl`
- `api` → [[gncz.station_geometry_rpostationgeometry|RPOStationGeometry]] · `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/station_geometry.jl`
- `api` → [[gncz.thruster_guidance_functions__osculating_elements_and_periapsis_direction|_osculating_elements_and_periapsis_direction]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl`
- `api` → [[gncz.thruster_guidance_models_aerobrakingcampaignpropulsivemaneuverguidancemodel|AerobrakingCampaignPropulsiveManeuverGuidanceModel]] · `module_api` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_models.jl`
- `api` → [[grp.src_gnc_command_types_jl|gnc/command_types.jl]] · `members_in` · call · `src/gnc/command_types.jl`
- `api` → [[grp.src_gnc_control|gnc/control/]] · `members_in` · call · `src/gnc/control/aerobraking/control_commands.jl`
- `api` → [[grp.src_gnc_guidance|gnc/guidance/]] · `members_in` · call · `src/gnc/guidance/aerobraking/common/heat_rate_models.jl`
- `api` → [[grp.src_gnc_hypr|gnc/hypr/]] · `members_in` · call · `src/gnc/hypr/hypr_utils.jl`
- `api` → [[grp.src_gnc_internal|gnc/internal/]] · `members_in` · call · `src/gnc/internal/bridge_helpers.jl`
- `api` → [[grp.src_gnc_navigation|gnc/navigation/]] · `members_in` · call · `src/gnc/navigation/rpo_nav/distances/mesh_distance.jl`
- `api` → [[grp.src_gnc_robotics|gnc/robotics/]] · `members_in` · call · `src/gnc/robotics/robot_arm_hypr/config.jl`
- `api` → [[module.core|SimulationModel]] · `gnc` · call · `src/core/simulation_model.jl:29-36`
<!-- vulcan:connections:end -->

## Limitations
The OSQP solve in `rpo_lqmpc_control` is configured with `eps_abs = 1.0e-4`, `eps_rel = 1.0e-4` and `max_iter = 1000`; when the status is neither `:Solved` nor `:Solved_inaccurate` the function silently returns a zero acceleration, so an infeasible box-constrained problem degrades to open-loop coasting with no diagnostic. `rpo_ref_preview` clamps the plan index to the last available column, so a plan shorter than the horizon freezes the reference at its final waypoint instead of extrapolating. `cloth_ik` throws an `ErrorException` only when the residual exceeds ten times `position_tol_m` after `max_iters`, meaning a solution between one and ten tolerances is returned as converged.

## Provenance
Mapped from `src/gnc/command_types.jl`, `src/gnc/control/control_hooks.jl`, `src/gnc/control/rpo_mpc/lqmpc.jl` and `src/gnc/control/rpo_mpc/rpo_mpc_control_model.jl`.
