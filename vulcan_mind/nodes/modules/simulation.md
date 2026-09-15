---
id: module.simulation
label: RuntimeServices
kind: module
source:
  file: src/simulation/runtime_services.jl
  symbol: RuntimeServices
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: Shared SPICE/GRAM lock plus SimulationEngine, SimulationCampaigns,
    callback, checkpoint, and runtime service entrypoints included by the package.
tags:
- module
charts:
- master
origin: agent
---

# RuntimeServices

## Purpose
`RuntimeServices` establishes the shared native-library synchronization primitive used by SpaceAGORA’s SPICE and GRAM integrations. The runtime source deliberately aliases `GRAM_LOCK` to `SPICE_LOCK`, because both libraries can expose the same CSPICE symbols and therefore share C-level mutable state. The broader simulation module family contains engine configuration, public simulation execution, callbacks, checkpoint handling, adaptive campaign routing, and Monte Carlo orchestration.

## Theory & Math
The lock is a concurrency invariant rather than a physical equation: every native call that touches the shared CSPICE state must be serialized through one critical section. Simulation execution itself integrates the configured ODE, while Monte Carlo execution samples a configured scenario repeatedly and aggregates the returned sample values.

## Model & Assumptions
Native SPICE and GRAM calls are assumed to be non-reentrant across the shared C implementation. The engine assumes its configuration supplies a valid problem, solver, callbacks, and output policy. Campaign code assumes sample functions are deterministic enough for the chosen randomization and that worker processes can reconstruct their configuration.

## Design & Implementation
`runtime_services.jl` defines `SPICE_LOCK` and the alias `GRAM_LOCK`. The simulation aggregator includes engine configuration and public API files, callbacks, checkpointing, campaign routing, and Monte Carlo support. `run_simulation` owns the single-scenario execution path; `run_monte_carlo` builds repeated samples from `MonteCarloSpec`, combines route features, and returns `MonteCarloResult`. The parallel profile namespace supplies optional worker routing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | Shared SPICE/GRAM lock plus SimulationEngine, SimulationCampaigns, callback, checkpoint, and runtime service entrypoints included by the package. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_in` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[grp.src_simulation_campaigns|simulation/campaigns/]] · `members_in` · call · `src/simulation/campaigns/adaptive_routing.jl`
- `api` → [[grp.src_simulation_engine|simulation/engine/]] · `members_in` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[grp.src_simulation_runtime_services_jl|simulation/runtime_services.jl]] · `members_in` · call · `src/simulation/runtime_services.jl`
- `api` → [[module.spaceagora|SpaceAGORA]] · `simulation` · call · `src/SpaceAGORA.jl:8-11`
- `api` → [[simulation.adaptive_routing__campaign_features_for_routing|_campaign_features_for_routing]] · `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- `api` → [[simulation.adaptive_routing__campaign_route_plan|_campaign_route_plan]] · `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- `api` → [[simulation.adaptive_routing__record_campaign_route_feedback_bang|_record_campaign_route_feedback!]] · `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- `api` → [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `module_api` · call · `src/simulation/campaigns/adaptive_routing.jl`
- `api` → [[simulation.assembly__callback_tolerances_for_phase|_callback_tolerances_for_phase]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[simulation.assembly__entry_target_count|_entry_target_count]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[simulation.assembly__requires_density_callback|_requires_density_callback]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[simulation.assembly__requires_density_for_rhs|_requires_density_for_rhs]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[simulation.assembly__requires_guidance_orbit_counter|_requires_guidance_orbit_counter]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[simulation.assembly__resolved_component_tolerance|_resolved_component_tolerance]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[simulation.assembly__uses_atmospheric_dynamic_effector|_uses_atmospheric_dynamic_effector]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`
- `api` → [[simulation.config__callback_outer_parallel_hint|_callback_outer_parallel_hint]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`
- `api` → [[simulation.config__control_callback_use_threads|_control_callback_use_threads]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`
- `api` → [[simulation.config__density_callback_use_threads|_density_callback_use_threads]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`
- `api` → [[simulation.config__policy_env_config|_policy_env_config]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`
- `api` → [[simulation.config_gramtrackcache|GramTrackCache]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/config.jl`
- `api` → [[simulation.constellation_ensemble__ensemble_member_settings|_ensemble_member_settings]] · `module_api` · call · `src/simulation/campaigns/constellation_ensemble.jl`
- `api` → [[simulation.control_callbacks__maybe_add_control_tstop_bang|_maybe_add_control_tstop!]] · `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- `api` → [[simulation.control_callbacks__run_guidance_for_thruster_schedule_bang|_run_guidance_for_thruster_schedule!]] · `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- `api` → [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- `api` → [[simulation.control_callbacks_schedule_all_bang|schedule_all!]] · `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- `api` → [[simulation.control_callbacks_schedule_event_driven_thruster_controls_bang|schedule_event_driven_thruster_controls!]] · `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_batchable_effector_flat_bang|_accumulate_batchable_effector_flat!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_flat_batch_bang|_accumulate_dynamic_effectors_flat_batch!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_dynamic_effectors_partitioned_bang|_accumulate_dynamic_effectors_partitioned!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_harmonics_flat_batch_bang|_accumulate_harmonics_flat_batch!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_invsq_flat_batch_bang|_accumulate_invsq_flat_batch!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_invsq_j2_flat_batch_bang|_accumulate_invsq_j2_flat_batch!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_nbody_flat_batch_bang|_accumulate_nbody_flat_batch!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__accumulate_srp_flat_batch_bang|_accumulate_srp_flat_batch!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__all_active_spacecraft_outside_atmosphere|_all_active_spacecraft_outside_atmosphere]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__any_robot_arm_effector|_any_robot_arm_effector]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__build_constellation_execution_plan_bang|_build_constellation_execution_plan!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__constellation_node_eff_idx|_constellation_node_eff_idx]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__constellation_node_sat_idx|_constellation_node_sat_idx]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__constellation_node_work_item|_constellation_node_work_item]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__count_flat_queue_only_effectors|_count_flat_queue_only_effectors]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__count_non_batchable_effectors|_count_non_batchable_effectors]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__drag_state_buffer_current|_drag_state_buffer_current]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__effector_value|_effector_value]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__ensure_rhs_flat_effector_scratch_bang|_ensure_rhs_flat_effector_scratch!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__flat_partition_selected|_flat_partition_selected]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__flat_totals_force_torque|_flat_totals_force_torque]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__gravity_backbone_core_acceleration|_gravity_backbone_core_acceleration]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__gravity_backbone_half_kick_bang|_gravity_backbone_half_kick!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__gravity_backbone_kick_acceleration|_gravity_backbone_kick_acceleration]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__gravity_backbone_state_sample|_gravity_backbone_state_sample]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__gravity_backbone_xyz_chunk|_gravity_backbone_xyz_chunk]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__has_any_batchable_effector|_has_any_batchable_effector]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__has_any_harmonics_effector|_has_any_harmonics_effector]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__partition_needs_state_sample|_partition_needs_state_sample]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__partition_selected_count|_partition_selected_count]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__prefill_atmosphere_samples_bang|_prefill_atmosphere_samples!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__prefill_rhs_flat_state_samples_bang|_prefill_rhs_flat_state_samples!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__prepare_rhs_flat_work_items_bang|_prepare_rhs_flat_work_items!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__prepare_rhs_flat_work_packets_bang|_prepare_rhs_flat_work_packets!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_effector_cost_rank|_rhs_effector_cost_rank]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_effector_estimated_cost_ns|_rhs_effector_estimated_cost_ns]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_effector_static_cost_ns|_rhs_effector_static_cost_ns]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_flat_item_eff_idx|_rhs_flat_item_eff_idx]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_flat_item_estimated_cost_ns|_rhs_flat_item_estimated_cost_ns]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_flat_packet_work_stats|_rhs_flat_packet_work_stats]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_flat_state_sample_from_buffers|_rhs_flat_state_sample_from_buffers]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__rhs_flat_use_packet_scheduler|_rhs_flat_use_packet_scheduler]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__robot_arm_coupling|_robot_arm_coupling]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__robot_arm_coupling_from_effector|_robot_arm_coupling_from_effector]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__robot_arm_effector_matches|_robot_arm_effector_matches]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__robot_arm_present|_robot_arm_present]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__solver_partition_validated|_solver_partition_validated]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__spacecraft_outside_atmosphere_for_current_state|_spacecraft_outside_atmosphere_for_current_state]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__update_rhs_flat_packet_cost_model_bang|_update_rhs_flat_packet_cost_model!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__update_rhs_flat_packet_overhead_model_bang|_update_rhs_flat_packet_overhead_model!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs__with_packet_scheduler|_with_packet_scheduler]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs_constellationexecutionplan|ConstellationExecutionPlan]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs_constellationinteractionedgeworkitem|ConstellationInteractionEdgeWorkItem]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs_spacecraft_dynamics_fast_control_bang|spacecraft_dynamics_fast_control!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs_spacecraft_dynamics_gravity_backbone_bang|spacecraft_dynamics_gravity_backbone!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.dynamics_rhs_spacecraft_dynamics_slow_bang|spacecraft_dynamics_slow!]] · `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`
- `api` → [[simulation.effector_sampling__buffered_atmosphere_valid|_buffered_atmosphere_valid]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling__extract_sample_mass_kg|_extract_sample_mass_kg]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling__extract_sample_pos_vel|_extract_sample_pos_vel]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling__planet_lpi_at_engine|_planet_lpi_at_engine]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling__sample_reusable_atmosphere|_sample_reusable_atmosphere]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling__sample_reusable_planet_frame|_sample_reusable_planet_frame]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling__sample_reusable_solar|_sample_reusable_solar]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling__wrench_method_available|_wrench_method_available]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling_build_state_sample|build_state_sample]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling_sample_atmosphere|sample_atmosphere]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling_sample_buffered_atmosphere|sample_buffered_atmosphere]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling_sample_buffered_planet_frame|sample_buffered_planet_frame]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling_sample_environment_with_buffered_atm|sample_environment_with_buffered_atm]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling_sample_environment_with_reusable_buffers|sample_environment_with_reusable_buffers]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.effector_sampling_sample_planet_frame_with_lpi|sample_planet_frame_with_lpi]] · `module_api` · call · `src/simulation/engine/effector_sampling.jl`
- `api` → [[simulation.event_callbacks_affect_upcrossing_bang|affect_upcrossing!]] · `module_api` · call · `src/simulation/callbacks/event_callbacks.jl`
- `api` → [[simulation.event_callbacks_save_func|save_func]] · `module_api` · call · `src/simulation/callbacks/event_callbacks.jl`
- `api` → [[simulation.from_env__engine_env_get_with_env_fallback|_engine_env_get_with_env_fallback]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__engine_env_haskey|_engine_env_haskey]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__engine_env_haskey_with_env_fallback|_engine_env_haskey_with_env_fallback]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__engine_env_overrides|_engine_env_overrides]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__env_bool|_env_bool]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__parse_float_opt|_parse_float_opt]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__parse_multirate_solver_sym|_parse_multirate_solver_sym]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__parse_solver_mode_sym|_parse_solver_mode_sym]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.from_env__parse_split_imex_solver_sym|_parse_split_imex_solver_sym]] · `module_api` · call · `src/simulation/engine/adapters/from_env.jl`
- `api` → [[simulation.interpolation__gram_track_cache_eval|_gram_track_cache_eval]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl`
- `api` → [[simulation.interpolation__gram_track_cache_profile|_gram_track_cache_profile]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl`
- `api` → [[simulation.model_selection__gram_density_cache_for_sat_bang|_gram_density_cache_for_sat!]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl`
- `api` → [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- `api` → [[simulation.monte_carlo__run_monte_carlo_sample|_run_monte_carlo_sample]] · `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- `api` → [[simulation.monte_carlo__throw_first_monte_carlo_failure|_throw_first_monte_carlo_failure]] · `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- `api` → [[simulation.monte_carlo_montecarlosampleresult|MonteCarloSampleResult]] · `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- `api` → [[simulation.persistence__append_saved_segment_bang|_append_saved_segment!]] · `module_api` · call · `src/simulation/engine/persistence.jl`
- `api` → [[simulation.persistence__append_series_columns_bang|_append_series_columns!]] · `module_api` · call · `src/simulation/engine/persistence.jl`
- `api` → [[simulation.persistence__collision_results_csv_path|_collision_results_csv_path]] · `module_api` · call · `src/simulation/engine/persistence.jl`
- `api` → [[simulation.persistence__find_sample_value|_find_sample_value]] · `module_api` · call · `src/simulation/engine/persistence.jl`
- `api` → [[simulation.persistence__results_bundle_prefix|_results_bundle_prefix]] · `module_api` · call · `src/simulation/engine/persistence.jl`
- `api` → [[simulation.persistence__results_csv_path|_results_csv_path]] · `module_api` · call · `src/simulation/engine/persistence.jl`
- `api` → [[simulation.planet_frame__planet_lpi_from_backend|_planet_lpi_from_backend]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl`
- `api` → [[simulation.planet_frame__planet_lpi_from_cache|_planet_lpi_from_cache]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl`
- `api` → [[simulation.planet_frame__planet_relative_state|_planet_relative_state]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl`
- `api` → [[simulation.public_api_load_nbody_ephemeris_cache_bang|load_nbody_ephemeris_cache!]] · `module_api` · call · `src/simulation/engine/public_api.jl`
- `api` → [[simulation.public_api_prewarm_nbody_ephemeris_cache|prewarm_nbody_ephemeris_cache]] · `module_api` · call · `src/simulation/engine/public_api.jl`
- `api` → [[simulation.registry__gram_runtime_stats_reset_bang|_gram_runtime_stats_reset!]] · `module_api` · call · `src/simulation/callbacks/registry.jl`
- `api` → [[simulation.registry__gram_runtime_stats_snapshot|_gram_runtime_stats_snapshot]] · `module_api` · call · `src/simulation/callbacks/registry.jl`
- `api` → [[simulation.resume_checkpoint__checkpoint_directory|_checkpoint_directory]] · `module_api` · call · `src/simulation/engine/resume_checkpoint.jl`
- `api` → [[simulation.resume_checkpoint__checkpoint_paths|_checkpoint_paths]] · `module_api` · call · `src/simulation/engine/resume_checkpoint.jl`
- `api` → [[simulation.resume_checkpoint__clear_checkpoint_bang|_clear_checkpoint!]] · `module_api` · call · `src/simulation/engine/resume_checkpoint.jl`
- `api` → [[simulation.rhs_calibration__calib_machine_label|_calib_machine_label]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__calib_sat_bucket|_calib_sat_bucket]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__make_calib_flat_plan|_make_calib_flat_plan]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__make_calib_satellite_batch_plan|_make_calib_satellite_batch_plan]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__rhs_calib_load_bang|_rhs_calib_load!]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__rhs_calib_path|_rhs_calib_path]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__rhs_calibrate_n_timed|_rhs_calibrate_n_timed]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__rhs_calibrate_n_warmup|_rhs_calibrate_n_warmup]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.rhs_calibration__rhs_plan_candidates|_rhs_plan_candidates]] · `module_api` · call · `src/simulation/engine/rhs_calibration.jl`
- `api` → [[simulation.runtime__buffered_stage_environment_state|_buffered_stage_environment_state]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`
- `api` → [[simulation.runtime__density_segment_end_t|_density_segment_end_t]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`
- `api` → [[simulation.runtime__extract_pos_vel|_extract_pos_vel]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`
- `api` → [[simulation.runtime__stage_environment_state|_stage_environment_state]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/runtime.jl`
- `api` → [[simulation.save_fields__save_snapshot|_save_snapshot]] · `module_api` · call · `src/simulation/callbacks/save_fields.jl`
- `api` → [[simulation.setup__any_effector_consumes_atmosphere|_any_effector_consumes_atmosphere]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__body_query_names_reuse_key|_body_query_names_reuse_key]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__cache_from_nbody_ephemeris_payload|_cache_from_nbody_ephemeris_payload]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__cache_time_key|_cache_time_key]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__clear_ephemeris_reuse_cache_bang|_clear_ephemeris_reuse_cache!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__collect_nbody_query_names|_collect_nbody_query_names]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__dynamic_effector_threadsafe|_dynamic_effector_threadsafe]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__dynamic_effectors_parallel_supported|_dynamic_effectors_parallel_supported]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_allow_with_outer|_effector_allow_with_outer]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_cost_ema_alpha|_effector_cost_ema_alpha]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_cost_min_samples|_effector_cost_min_samples]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_cost_ns_per_item_default|_effector_cost_ns_per_item_default]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_heavy_only|_effector_heavy_only]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_long_mission_threshold_s|_effector_long_mission_threshold_s]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_long_orbit_threshold|_effector_long_orbit_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_max_threads|_effector_max_threads]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_observed_cost_ns_per_item|_effector_observed_cost_ns_per_item]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_outer_parallel_hint|_effector_outer_parallel_hint]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_outer_work_scale|_effector_outer_work_scale]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_parallel_mode|_effector_parallel_mode]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_satellite_share_budget|_effector_satellite_share_budget]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_shared_buffers|_effector_shared_buffers]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_thread_threshold|_effector_thread_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__effector_work_ns_per_worker_threshold|_effector_work_ns_per_worker_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ensure_rhs_effector_cost_model_bang|_ensure_rhs_effector_cost_model!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ephemerides_model_reuse_key|_ephemerides_model_reuse_key]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ephemerides_time_seconds_flexible|_ephemerides_time_seconds_flexible]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ephemeris_explicit_cache_store_bang|_ephemeris_explicit_cache_store!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ephemeris_reuse_enabled|_ephemeris_reuse_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ephemeris_reuse_lookup|_ephemeris_reuse_lookup]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ephemeris_reuse_max_entries|_ephemeris_reuse_max_entries]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__ephemeris_reuse_store_bang|_ephemeris_reuse_store!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__gram_per_sat_instances_enabled|_gram_per_sat_instances_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__harmonics_batch_allow_with_outer|_harmonics_batch_allow_with_outer]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__has_active_nbody_effector|_has_active_nbody_effector]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__has_active_srp_effector|_has_active_srp_effector]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__is_nbody_effector_like|_is_nbody_effector_like]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__load_nbody_ephemeris_cache_bang|_load_nbody_ephemeris_cache!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__mission_is_long_for_effector_threads|_mission_is_long_for_effector_threads]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__nbody_ephemeris_body_index_by_name|_nbody_ephemeris_body_index_by_name]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__nbody_ephemeris_cache_dt_s|_nbody_ephemeris_cache_dt_s]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__nbody_ephemeris_cache_enabled|_nbody_ephemeris_cache_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__nbody_ephemeris_cache_max_samples|_nbody_ephemeris_cache_max_samples]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__nbody_ephemeris_cache_payload|_nbody_ephemeris_cache_payload]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__nbody_ephemeris_reuse_key|_nbody_ephemeris_reuse_key]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__parse_nonnegative_int_env|_parse_nonnegative_int_env]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__parse_unit_float_env|_parse_unit_float_env]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__payload_field|_payload_field]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__planet_frame_cache_dt_s|_planet_frame_cache_dt_s]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__planet_frame_cache_enabled|_planet_frame_cache_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__planet_frame_cache_max_samples|_planet_frame_cache_max_samples]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__planet_frame_ephemeris_reuse_key|_planet_frame_ephemeris_reuse_key]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__planet_transform_key|_planet_transform_key]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__policy_env_config|_policy_env_config]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__prewarmed_nbody_ephemeris_lookup|_prewarmed_nbody_ephemeris_lookup]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__profile_forces_serial_rhs|_profile_forces_serial_rhs]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__register_prewarmed_nbody_ephemeris_cache_bang|_register_prewarmed_nbody_ephemeris_cache!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_batch_parallel_mode|_rhs_batch_parallel_mode]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_batch_thread_threshold|_rhs_batch_thread_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_effector_cost_min_samples|_rhs_effector_cost_min_samples]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_effector_observed_cost_ns|_rhs_effector_observed_cost_ns]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_effectors_have_heavy_or_heterogeneous_cost|_rhs_effectors_have_heavy_or_heterogeneous_cost]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_execution_mode_env|_rhs_execution_mode_env]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_batch_privileged_effector|_rhs_flat_batch_privileged_effector]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_cost_heterogeneity_threshold|_rhs_flat_cost_heterogeneity_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_has_batch_privileged_effector|_rhs_flat_has_batch_privileged_effector]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_min_effectors|_rhs_flat_min_effectors]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_min_sats|_rhs_flat_min_sats]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_min_thread_budget|_rhs_flat_min_thread_budget]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_packet_heterogeneity_threshold|_rhs_flat_packet_heterogeneity_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_packet_min_items|_rhs_flat_packet_min_items]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_packet_overhead_disable_ratio|_rhs_flat_packet_overhead_disable_ratio]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_packet_overhead_min_samples|_rhs_flat_packet_overhead_min_samples]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_packet_scheduler_mode|_rhs_flat_packet_scheduler_mode]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_packet_target_min_ns|_rhs_flat_packet_target_min_ns]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_packet_work_ns_threshold|_rhs_flat_packet_work_ns_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_supported|_rhs_flat_supported]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_work_ns_threshold|_rhs_flat_work_ns_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_flat_work_per_worker_ns_threshold|_rhs_flat_work_per_worker_ns_threshold]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_harmonics_batch_enabled|_rhs_harmonics_batch_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_harmonics_batch_min_sats_per_worker|_rhs_harmonics_batch_min_sats_per_worker]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_harmonics_flat_experimental_enabled|_rhs_harmonics_flat_experimental_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_invsq_flat_min_sats|_rhs_invsq_flat_min_sats]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_plan_step_cache_enabled|_rhs_plan_step_cache_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_single_harmonics_flat_supported|_rhs_single_harmonics_flat_supported]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__rhs_single_invsq_flat_supported|_rhs_single_invsq_flat_supported]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__satellite_batch_saturates_pool|_satellite_batch_saturates_pool]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__snapshot_rhs_plan_env_config|_snapshot_rhs_plan_env_config]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__spice_rhs_memo_enabled|_spice_rhs_memo_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__srp_ephemeris_cache_dt_s|_srp_ephemeris_cache_dt_s]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__srp_ephemeris_cache_enabled|_srp_ephemeris_cache_enabled]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__srp_ephemeris_cache_max_samples|_srp_ephemeris_cache_max_samples]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__srp_ephemeris_reuse_key|_srp_ephemeris_reuse_key]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__update_effector_cost_model_bang|_update_effector_cost_model!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__update_rhs_effector_cost_model_bang|_update_rhs_effector_cost_model!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__with_serial_effector_decision|_with_serial_effector_decision]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.setup__write_nbody_ephemeris_cache_file_bang|_write_nbody_ephemeris_cache_file!]] · `module_api` · call · `src/simulation/engine/setup.jl`
- `api` → [[simulation.solver_policy__auto_stiff_smooth_gravity_effector|_auto_stiff_smooth_gravity_effector]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__auto_stiff_smooth_gravity_reject_reason|_auto_stiff_smooth_gravity_reject_reason]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__cache_integrator_bang|_cache_integrator!]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__gravity_backbone_fixed_dt_s|_gravity_backbone_fixed_dt_s]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__gravity_backbone_has_kicks|_gravity_backbone_has_kicks]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__gravity_backbone_kick_structure_validated|_gravity_backbone_kick_structure_validated]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__gravity_backbone_structure_validated|_gravity_backbone_structure_validated]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__multirate_slow_dt_s|_multirate_slow_dt_s]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__requires_componentwise_tolerances|_requires_componentwise_tolerances]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__resolve_component_tolerance|_resolve_component_tolerance]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__retcode_is_stiff_symptom|_retcode_is_stiff_symptom]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__solver_bool_env|_solver_bool_env]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__solver_cache_options_match|_solver_cache_options_match]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__solver_save_everystep|_solver_save_everystep]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy__split_subproblem|_split_subproblem]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.solver_policy_solverintegratorcache|SolverIntegratorCache]] · `module_api` · call · `src/simulation/engine/solver_policy.jl`
- `api` → [[simulation.state_access__gravity_backbone_initial_states|_gravity_backbone_initial_states]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__gravity_backbone_position_state|_gravity_backbone_position_state]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__gravity_backbone_velocity_state|_gravity_backbone_velocity_state]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__state_has_heat_loads|_state_has_heat_loads]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__state_has_mass|_state_has_mass]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__state_has_quaternion|_state_has_quaternion]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__state_heat_loads|_state_heat_loads]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__state_mass_kg|_state_mass_kg]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__state_quaternion|_state_quaternion]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.state_access__state_velocity_ii|_state_velocity_ii]] · `module_api` · call · `src/simulation/engine/state_access.jl`
- `api` → [[simulation.targeting__eccentric_to_true_anomaly|_eccentric_to_true_anomaly]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- `api` → [[simulation.targeting__gram_kepler_target|_gram_kepler_target]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- `api` → [[simulation.targeting__gram_linear_target|_gram_linear_target]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- `api` → [[simulation.targeting__solve_kepler_elliptic|_solve_kepler_elliptic]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- `api` → [[simulation.targeting__true_to_eccentric_anomaly|_true_to_eccentric_anomaly]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl`
- `api` → [[simulation.thermal_callbacks__heat_rate_buffer_for_sat_bang|_heat_rate_buffer_for_sat!]] · `module_api` · call · `src/simulation/callbacks/thermal_callbacks.jl`
- `api` → [[simulation.vacuum_predicted_gram__eval_natural_cubic_spline|_eval_natural_cubic_spline]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation.vacuum_predicted_gram__interp_vacuum_alt|_interp_vacuum_alt]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation.vacuum_predicted_gram__interp_vacuum_position|_interp_vacuum_position]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation.vacuum_predicted_gram__interp_vacuum_wind|_interp_vacuum_wind]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation.vacuum_predicted_gram__vacuum_gram_cache_for_sat_bang|_vacuum_gram_cache_for_sat!]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation.vacuum_predicted_gram__vacuum_j2_accel|_vacuum_j2_accel]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation.vacuum_predicted_gram_vacuumpredictedgramcache|VacuumPredictedGRAMCache]] · `module_api` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl`
- `api` → [[simulation_a.callbacks_simulationcallbacks|SimulationCallbacks]] · `module_api` · call · `src/simulation/callbacks/callbacks.jl`
- `api` → [[simulation_a.density_callbacks_density_callbacks|density_callbacks]] · `module_api` · call · `src/simulation/callbacks/density_callbacks.jl`
- `api` → [[simulation_a.gram_track_cache_gram_track_cache|gram_track_cache]] · `module_api` · call · `src/simulation/callbacks/gram_track_cache.jl`
- `api` → [[simulation_a.registry_gramruntimestats|GramRuntimeStats]] · `module_api` · call · `src/simulation/callbacks/registry.jl`
- `api` → [[simx.campaigns_simulation_campaigns_simulationcampaigns|SimulationCampaigns]] · `module_api` · call · `src/simulation/campaigns/simulation_campaigns.jl`
- `api` → [[simx.engine_persistence_write_results_bundle__write_results_bundle_bang|_write_results_bundle!]] · `module_api` · call · `src/simulation/engine/persistence.jl`
- `api` → [[simx.engine_simulation_engine_simulationengine|SimulationEngine]] · `module_api` · call · `src/simulation/engine/simulation_engine.jl`
- `api` → [[simx.runtime_services_runtimeservices|RuntimeServices]] · `module_api` · call · `src/simulation/runtime_services.jl`
<!-- vulcan:connections:end -->

## Limitations
The shared lock prevents the documented native race but does not make arbitrary external library calls safe when they bypass SpaceAGORA’s lock. Solver failure, invalid configuration, or a non-serializable sample closure can terminate a run. Monte Carlo output quality depends on the scenario generator and number of samples; the runtime does not infer statistical convergence automatically.

## Provenance
Mapped from `src/simulation/runtime_services.jl` and the included files under `src/simulation/engine`, `callbacks`, `campaigns`, and `checkpoint`.
