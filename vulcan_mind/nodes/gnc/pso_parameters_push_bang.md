---
id: gnc.pso_parameters_push_bang
label: push!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: push!
  lines:
  - 440
  - 440
inputs:
- id: pairs
  type: Any
  units: n/a
  required: true
  description: Positional argument `pairs`.
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
  description: Return value of `push!`; mutates `pairs` in place. Returns `> value)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# push!

## Purpose
This node records the `push!` call site at line 440 of `pso_parameters.jl`, inside `_rpo_pso_normalize_kwargs`, where each incoming keyword pair is appended to the `pairs` vector after alias translation. It is `Base.push!`, not a method defined by this file.

## Design & Implementation
The expression is `push!(pairs, get(RPO_PSO_CONFIG_ALIASES, key, key) => value)` with `pairs::Vector{Pair{Symbol, Any}}`. `push!` grows the vector in place by one element per keyword, preserving caller order so that later duplicates override earlier ones when the vector is splatted into a NamedTuple via `(; pairs...)`. The mutated object is only the local `pairs` vector, which never escapes the function.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pairs` | Any | n/a | yes | Positional argument `pairs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `push!`; mutates `pairs` in place. Returns `> value)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.calibration__annotate_calibration_rows|_annotate_calibration_rows]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:84-84`
- [[analysis.comparison_metrics__rates|_rates]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl:83-83`
- [[analysis.decay_diagnostics_flight_density_table|flight_density_table]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:160-160`
- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:580-580`
- [[analysis.manifest_parsing__optional_symbol_vector|_optional_symbol_vector]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:103-103`
- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:383-383`
- [[analysis.scenario_builders__make_tabulated_flight_density_model|_make_tabulated_flight_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:221-221`
- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:131-131`
- [[analysis.scenario_builders__telemetry_coefficients_normalized_for_scenario|_telemetry_coefficients_normalized_for_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:83-83`
- [[analysis.scenario_builders__telemetry_j2_source_for_scenario|_telemetry_j2_source_for_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:68-68`
- [[analysis.telemetry_loading__extract_extrema_from_time_aligned_telemetry|_extract_extrema_from_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:464-464`
- [[analysis.telemetry_loading__extract_extrema_series|_extract_extrema_series]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:128-128`
- [[assets.rpo_station_assets__load_stl_triangles|_load_stl_triangles]] · `callees` → `callers` · call · `src/assets/rpo_station_assets.jl:88-88`
- [[cli.assets_load_asset_manifest|load_asset_manifest]] · `callees` → `callers` · call · `src/cli/assets.jl:52-52`
- [[cli.spaceagora_cli__run_benchmark|_run_benchmark]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:160-160`
- [[cli.spaceagora_cli__run_example|_run_example]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:94-94`
- [[dynamics.cloth_multibody_build_compliant_topology|build_compliant_topology]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:303-303`
- [[dynamics.cloth_multibody_build_rectangular_compliant_grid|build_rectangular_compliant_grid]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:375-375`
- [[dynamics.cloth_multibody_compliant_joint_loads|compliant_joint_loads]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:562-562`
- [[dynamics.cloth_multibody_simulate_compliant_multibody|simulate_compliant_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:715-715`
- [[dynamics.cloth_robot_arm_dynamics_cloth_robot_arm_multibody|cloth_robot_arm_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:257-257`
- [[dynamics.cloth_robot_arm_dynamics_simulate_cloth_robot_arm_plan|simulate_cloth_robot_arm_plan]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_robot_arm_dynamics.jl:470-470`
- [[dynamics.perturbations__harmonics_model_cache_key|_harmonics_model_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:886-886`
- [[envana.ana_decay_diagnostics_secular_sma_slope|secular_sma_slope]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:61-61`
- [[envana.ana_rpo_visualization_rpo_path_plot|rpo_path_plot]] · `callees` → `callers` · call · `src/analysis/visualization/rpo/rpo_visualization.jl:10-10`
- [[environment.planets__furnsh_once|_furnsh_once]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:47-47`
- [[gnc.constraint_tracking_time_switch_func_affect_bang|time_switch_func_affect!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:462-462`
- [[gnc.eom_predictor_shooting_residual_bang|shooting_residual!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1072-1072`
- [[gnc.hypr_utils_hypr_rrt_near_indices|hypr_rrt_near_indices]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:160-160`
- [[gnc.hypr_utils_hypr_rrt_refresh_subtree_costs_bang|hypr_rrt_refresh_subtree_costs!]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:205-205`
- [[gnc.hypr_utils_hypr_rrt_tree_path|hypr_rrt_tree_path]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:181-181`
- [[gnc.lqmpc_rpo_prediction_mats|rpo_prediction_mats]] · `callees` → `callers` · call · `src/gnc/control/rpo_mpc/lqmpc.jl:42-42`
- [[gnc.path_retiming_rpo_remove_near_duplicate_samples|rpo_remove_near_duplicate_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:26-26`
- [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:128-128`
- [[gnc.path_sampling_rpo_sample_path_polyline_adaptive|rpo_sample_path_polyline_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:181-181`
- [[gnc.planner_comparison_rpo_comparison_cost_iteration_plot|rpo_comparison_cost_iteration_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1084-1084`
- [[gnc.planner_comparison_rpo_comparison_failed_paths_plot|rpo_comparison_failed_paths_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1000-1000`
- [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:764-764`
- [[gnc.planner_comparison_rpo_comparison_path_family_plot|rpo_comparison_path_family_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:873-873`
- [[gnc.planner_comparison_rpo_comparison_single_path_plot|rpo_comparison_single_path_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:931-931`
- [[gnc.planner_comparison_rpo_group_metric_mean|rpo_group_metric_mean]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:683-683`
- [[gnc.planner_comparison_rpo_write_namedtuple_csv|rpo_write_namedtuple_csv]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:668-668`
- [[gnc.planner_core_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:289-289`
- [[gnc.pso_path_planning_record_iteration_bang|record_iteration!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:266-266`
- [[gnc.pso_path_planning_record_iteration_timeout_bang|record_iteration_timeout!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:276-276`
- [[gnc.robot_arm_planning__reference_times|_reference_times]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_planning.jl:64-64`
- [[gnc.rpo_guidance_hooks__rpo_record_replanning_event_bang|_rpo_record_replanning_event!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:68-68`
- [[gnc.rrt_connect_rpo_rrt_extend_bang|rpo_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:173-173`
- [[gnc.rrt_connect_rpo_rrt_star_add_node_bang|rpo_rrt_star_add_node!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:252-252`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:574-574`
- [[gnc.rrt_warmstart__robot_arm_rrt_extend_bang|_robot_arm_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:74-74`
- [[gnc.targeting_control__edg_certify_targeting_candidates|_edg_certify_targeting_candidates]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:909-909`
- [[gnc.thruster_guidance_functions__ensure_apo_target_state_bang|_ensure_apo_target_state!]] · `callees` → `callers` · call · `src/gnc/guidance/thruster_guidance/thruster_guidance_functions.jl:9-9`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:488-488`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:462-462`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1072-1072`
- [[gncy.path_retiming_rpo_retime_path|rpo_retime_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:193-193`
- [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:214-214`
- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:624-624`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:266-266`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:302-302`
- [[gncz.planner_core_plan_robot_arm_motion_hypr|plan_robot_arm_motion_hypr]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/planner_core.jl:289-289`
- [[grp.src_analysis_verification|analysis/verification/]] · `members_out` → `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:84-84`
- [[grp.src_analysis_visualization|analysis/visualization/]] · `members_out` → `callers` · call · `src/analysis/visualization/rpo/rpo_visualization.jl:10-10`
- [[grp.src_parallel_policy|parallel/policy/]] · `members_out` → `callers` · call · `src/parallel/policy/context.jl:302-302`
- [[grp.src_parallel_routing|parallel/routing/]] · `members_out` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:409-409`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:148-148`
- [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_out` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:197-197`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/execution.jl:78-78`
- [[grp.src_vehicle_robotics|vehicle/robotics/]] · `members_out` → `callers` · call · `src/vehicle/robotics/robotics.jl:189-189`
- [[grp.src_vehicle_spacecraft|vehicle/spacecraft/]] · `members_out` → `callers` · call · `src/vehicle/spacecraft/assembly.jl:59-59`
- [[grp.src_vehicle_structure|vehicle/structure/]] · `members_out` → `callers` · call · `src/vehicle/structure/mass_properties.jl:105-105`
- [[misc.assets_check_assets|check_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:88-88`
- [[misc.maneuver_plans_odyssey_campaign_maneuvers|odyssey_campaign_maneuvers]] · `callees` → `callers` · call · `src/mission/operations/maneuver_plans.jl:40-40`
- [[parallel.context__destroy_persistent_foreach_scope_bang|_destroy_persistent_foreach_scope!]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:302-302`
- [[parallel.outer_route_selection__route_ranked_candidates|_route_ranked_candidates]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:409-409`
- [[parallel.outer_route_selection_outer_route_candidates|outer_route_candidates]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:387-387`
- [[parallel.outer_route_state_load_outer_route_state_bang|load_outer_route_state!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:234-234`
- [[parallel.outer_route_state_save_outer_route_state|save_outer_route_state]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:163-163`
- [[parallel.persistent_hints__hint_candidate_allotments|_hint_candidate_allotments]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:214-214`
- [[parallel.persistent_hints__save_persistent_hint_state_locked_bang|_save_persistent_hint_state_locked!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:123-123`
- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:419-419`
- [[parcore.persistent_hints__hint_choose_allotment|_hint_choose_allotment]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:316-316`
- [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:197-197`
- [[simulation.constellation_ensemble__validate_ensemble_uncoupled|_validate_ensemble_uncoupled]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:55-55`
- [[simulation.execution__append_backbone_saved_segment_bang|_append_backbone_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:78-78`
- [[simulation.execution__append_checkpoint_saved_segment_bang|_append_checkpoint_saved_segment!]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:98-98`
- [[simulation.execution__build_block_diagonal_jac_prototype|_build_block_diagonal_jac_prototype]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:16-16`
- [[simulation.model_selection__ensure_gram_isolated_pool_bang|_ensure_gram_isolated_pool!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:148-148`
- [[simulation.monte_carlo__run_monte_carlo_serial|_run_monte_carlo_serial]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:113-113`
- [[simulation.rhs_calibration__rhs_calib_save_bang|_rhs_calib_save!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:135-135`
- [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:231-231`
- [[simulation.setup__collect_nbody_query_names|_collect_nbody_query_names]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1483-1483`
- [[simulation.setup__initialize_density_model_instances_bang|_initialize_density_model_instances!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1329-1329`
- [[simulation.solver_policy__solve_with_gravity_backbone_solver|_solve_with_gravity_backbone_solver]] · `callees` → `callers` · call · `src/simulation/engine/solver_policy.jl:592-592`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:128-128`
- [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callees` → `callers` · call · `src/simulation/callbacks/save_fields.jl:188-188`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:357-357`
- [[vehicle.assembly_add_facet_bang|add_facet!]] · `callees` → `callers` · call · `src/vehicle/spacecraft/assembly.jl:59-59`
- [[vehicle.assembly_add_joint_bang|add_joint!]] · `callees` → `callers` · call · `src/vehicle/spacecraft/assembly.jl:54-54`
- [[vehicle.assembly_add_magnet_bang|add_magnet!]] · `callees` → `callers` · call · `src/vehicle/spacecraft/assembly.jl:67-67`
- [[vehicle.assembly_add_thruster_bang|add_thruster!]] · `callees` → `callers` · call · `src/vehicle/spacecraft/assembly.jl:90-90`
- [[vehicle.cloth_fk|cloth_fk]] · `callees` → `callers` · call · `src/vehicle/robotics/robotics.jl:189-189`
- [[vehicle.mass_properties_get_spacecraft_mass|get_spacecraft_mass]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:105-105`
- [[vehx.spacecraft_assembly_add_body_bang|add_body!]] · `callees` → `callers` · call · `src/vehicle/spacecraft/assembly.jl:34-34`
- [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callees` → `callers` · call · `src/vehicle/structure/assembly_graph.jl:12-12`
- [[vehx.structure_geometry_properties_get_spacecraft_reference_area|get_spacecraft_reference_area]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:10-10`
- [[vehx.structure_mass_properties_update_inertia_tensor_bang|update_inertia_tensor!]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:41-41`

**Downstream**

- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:579-579`
- `callees` → [[gncy.pso_parameters_rpopsoconfig|RPOPSOConfig]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:446-446`
<!-- vulcan:connections:end -->

## Limitations
Because this is a standard-library call rather than a local definition, there is no bespoke error handling; `push!` may reallocate the backing array as it grows, which is negligible for the tens of keywords typical here. The chart lists it as a symbol only because the mapper detected the call, so no separate interface exists beyond `Base.push!`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 440.
