---
id: module.analysis
label: TelemetryVerification
kind: module
source:
  file: src/analysis/verification/telemetry_verification.jl
  symbol: TelemetryVerification
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: VerificationRequest, VerificationResult, run_verification, run_verification_cli,
    and run_study exported by the telemetry verification package.
tags:
- module
charts:
- master
origin: agent
---

# TelemetryVerification

## Purpose
`TelemetryVerification` is the analysis boundary for comparing SpaceAGORA simulation output with telemetry-oriented benchmark scenarios. Its module declaration exports the request and result records together with the batch entrypoints used by command-line and study workflows. The source fixes repository-relative defaults for the output directory and benchmark manifest, then includes typed request definitions, manifest parsing, scenario construction, telemetry loading, comparison metrics, decay diagnostics, calibration, error tables, reporting, the runner, and inclination/covariance fitting. This ordering makes the runner depend on concrete helpers while keeping the public package surface small.

## Theory & Math
The verification tolerances are expressed as relative and absolute thresholds for orbital and atmospheric quantities. A comparison accepts a value when its residual is bounded by `atol + rtol * abs(reference)`. The module also records fixed timestep policies: 60 seconds for orbital telemetry and 0.2 seconds for atmospheric telemetry. These constants are part of the verification contract rather than solver parameters.

## Model & Assumptions
The default manifest is under `test/telemetry_benchmark_manifest.toml`, and default reports are written below the repository `output` directory. Verification runs through `SimulationEngine`, so the analysis layer assumes the same model construction and numerical conventions as the production simulation path. The strict tolerances are intended for regression comparison, not for identifying physical uncertainty.

## Design & Implementation
`telemetry_verification.jl` imports `SimulationModel` and `SimulationEngine`, includes the reference-system contract, and then assembles thirteen implementation files. `run_verification` and `run_study` are defined in the included runner and reporting layers, while manifest and scenario helpers isolate file-format concerns from metric calculation. `run_verification_cli` provides the command-facing adapter without changing the underlying request/result types.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | VerificationRequest, VerificationResult, run_verification, run_verification_cli, and run_study exported by the telemetry verification package. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[analysis.calibration__annotate_calibration_rows|_annotate_calibration_rows]] · `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`
- `api` → [[analysis.calibration__calibration_active|_calibration_active]] · `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`
- `api` → [[analysis.calibration__calibration_score|_calibration_score]] · `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`
- `api` → [[analysis.calibration__grid_values|_grid_values]] · `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`
- `api` → [[analysis.calibration__single_point_calibration|_single_point_calibration]] · `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`
- `api` → [[analysis.comparison_metrics__compare_time_series|_compare_time_series]] · `module_api` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl`
- `api` → [[analysis.comparison_metrics__rates|_rates]] · `module_api` · call · `src/analysis/verification/telemetry_verification/comparison_metrics.jl`
- `api` → [[analysis.decay_diagnostics_flight_density_table|flight_density_table]] · `module_api` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl`
- `api` → [[analysis.decay_diagnostics_visviva_sma|visviva_sma]] · `module_api` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl`
- `api` → [[analysis.decay_diagnostics_zero_referenced_decay|zero_referenced_decay]] · `module_api` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl`
- `api` → [[analysis.error_tables__telemetry_altitude_km|_telemetry_altitude_km]] · `module_api` · call · `src/analysis/verification/telemetry_verification/error_tables.jl`
- `api` → [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `module_api` · call · `src/analysis/verification/telemetry_verification/error_tables.jl`
- `api` → [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- `api` → [[analysis.example_support__example_smoke_enabled|_example_smoke_enabled]] · `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- `api` → [[analysis.example_support__example_smoke_mission_time|_example_smoke_mission_time]] · `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- `api` → [[analysis.example_support__example_smoke_results_enabled|_example_smoke_results_enabled]] · `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- `api` → [[analysis.example_support__link_q|_link_q]] · `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- `api` → [[analysis.example_support_make_three_body_spacecraft|make_three_body_spacecraft]] · `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- `api` → [[analysis.example_support_run_and_report|run_and_report]] · `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`
- `api` → [[analysis.ic_fit__ic_fit_series|_ic_fit_series]] · `module_api` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl`
- `api` → [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_bool|_optional_bool]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_float64_vector|_optional_float64_vector]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_float|_optional_float]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_float_tuple|_optional_float_tuple]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_int64_vector|_optional_int64_vector]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_int|_optional_int]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_str|_optional_str]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_str_vector|_optional_str_vector]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__optional_symbol_vector|_optional_symbol_vector]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_attitude_q|_parse_attitude_q]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_calibration_config|_parse_calibration_config]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_element_frame|_parse_element_frame]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_event_tolerance|_parse_event_tolerance]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_gravity_model|_parse_gravity_model]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_ic_offset|_parse_ic_offset]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_initial_time|_parse_initial_time]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_maneuver_config|_parse_maneuver_config]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_orbit_altitude_mode|_parse_orbit_altitude_mode]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_positive_int_env|_parse_positive_int_env]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_reference_frame|_parse_reference_frame]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_spacecraft_config|_parse_spacecraft_config]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_time_aligned_comparison_mode|_parse_time_aligned_comparison_mode]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_tolerances|_parse_tolerances]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_truth_mask|_parse_truth_mask]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_units|_parse_units]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__parse_vec3|_parse_vec3]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__request_from_study_config|_request_from_study_config]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__require_float|_require_float]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__require_int|_require_int]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__require_key|_require_key]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__require_str|_require_str]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__require_str_vector|_require_str_vector]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__require_table|_require_table]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__resolve_repo_path|_resolve_repo_path]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__telemetry_solver_maxiters|_telemetry_solver_maxiters]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__telemetry_solver_mode|_telemetry_solver_mode]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.manifest_parsing__telemetry_solver_retry_maxiters|_telemetry_solver_retry_maxiters]] · `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- `api` → [[analysis.reporting__add_scaled_column_bang|_add_scaled_column!]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__append_display_error_columns_bang|_append_display_error_columns!]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__append_display_metric_columns_bang|_append_display_metric_columns!]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__axis_units|_axis_units]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__default_plots_outdir|_default_plots_outdir]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__display_value_scale|_display_value_scale]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__display_value_units|_display_value_units]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__generate_plots|_generate_plots]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__maneuver_count|_maneuver_count]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__maneuver_replay_scale_mode|_maneuver_replay_scale_mode]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__orbit_altitude_mode|_orbit_altitude_mode]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__scenario_status_extra|_scenario_status_extra]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__source_file|_source_file]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.reporting__value_units|_value_units]] · `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`
- `api` → [[analysis.rpo_visualization_rpo_tracking_plot|rpo_tracking_plot]] · `module_api` · call · `src/analysis/visualization/rpo/rpo_visualization.jl`
- `api` → [[analysis.rpo_visualization_rpovisualization|RPOVisualization]] · `module_api` · call · `src/analysis/visualization/rpo/rpo_visualization.jl`
- `api` → [[analysis.runner__final_run_or_reused_eval|_final_run_or_reused_eval]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.runner__initial_condition_from_time_aligned_telemetry|_initial_condition_from_time_aligned_telemetry]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.runner__run_once|_run_once]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.runner__run_simulation_dataframe|_run_simulation_dataframe]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.runner__run_single_scenario|_run_single_scenario]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.runner__select_scenarios|_select_scenarios]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.runner_run_study|run_study]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.runner_run_verification_cli|run_verification_cli]] · `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`
- `api` → [[analysis.scenario_builders__base_gravity_effector|_base_gravity_effector]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__gramofflinesurrogatefallbackbase|_GRAMOfflineSurrogateFallbackBase]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__harmonics_order|_harmonics_order]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__has_campaign_maneuvers|_has_campaign_maneuvers]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__has_high_fidelity_effectors|_has_high_fidelity_effectors]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__is_gram_library_missing_error|_is_gram_library_missing_error]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__libraryless_gram_surrogate_enabled|_libraryless_gram_surrogate_enabled]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__make_required_gram_density_model|_make_required_gram_density_model]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__make_spacecraft|_make_spacecraft]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__make_tabulated_flight_density_model|_make_tabulated_flight_density_model]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__make_time_tabulated_density_model|_make_time_tabulated_density_model]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__nbody_primary_name|_nbody_primary_name]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__period_seconds|_period_seconds]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__planet_from_name|_planet_from_name]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__save_fields_for_study|_save_fields_for_study]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__scenario_initial_condition|_scenario_initial_condition]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__telemetry_coefficients_normalized_for_scenario|_telemetry_coefficients_normalized_for_scenario]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__telemetry_j2_source_for_scenario|_telemetry_j2_source_for_scenario]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__try_libraryless_gram_surrogate|_try_libraryless_gram_surrogate]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__with_environment_wind|_with_environment_wind]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__with_orbit_mission|_with_orbit_mission]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.scenario_builders_scaledaerodynamiccoefficientfm|ScaledAerodynamicCoefficientfM]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[analysis.telemetry_loading__differentiate_series|_differentiate_series]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__extract_extrema_from_time_aligned_telemetry|_extract_extrema_from_time_aligned_telemetry]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__initial_time_et|_initial_time_et]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__j2000_to_planet_fixed_state|_j2000_to_planet_fixed_state]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__load_time_aligned_telemetry|_load_time_aligned_telemetry]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__planet_fixed_frame_fallback_name|_planet_fixed_frame_fallback_name]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__planet_fixed_frame_name|_planet_fixed_frame_name]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__planet_fixed_to_j2000_state|_planet_fixed_to_j2000_state]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__require_column|_require_column]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.telemetry_loading__sun_unit_vector_j2000|_sun_unit_vector_j2000]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[analysis.types_abstractscenarioconfig|AbstractScenarioConfig]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[analysis.types_atmospheretruthconfig|AtmosphereTruthConfig]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[analysis.types_calibrationconfig|CalibrationConfig]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[analysis.types_orbiteventsscenarioconfig|OrbitEventsScenarioConfig]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[analysis.types_spacecraftconfig|SpacecraftConfig]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[analysis.types_timealignedscenarioconfig|TimeAlignedScenarioConfig]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[analysis.types_verificationrequest|VerificationRequest]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[analysis.types_verificationresult|VerificationResult]] · `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`
- `api` → [[envana.ana_scenario_builders_body_equator_frame_rotation|_body_equator_frame_rotation]] · `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`
- `api` → [[envana.ana_telemetry_loading_transform_state|_transform_state]] · `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`
- `api` → [[envana.ana_telemetry_verification_telemetryverification|TelemetryVerification]] · `module_api` · call · `src/analysis/verification/telemetry_verification.jl`
- `api` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `src/analysis/verification/telemetry_verification/calibration.jl`
- `api` → [[grp.src_analysis_visualization|analysis/visualization/]] · `members_in` · call · `src/analysis/visualization/rpo/rpo_visualization.jl`
- `api` → [[module.spaceagora|SpaceAGORA]] · `analysis` · call · `src/SpaceAGORA.jl:12-14`
<!-- vulcan:connections:end -->

## Limitations
The module depends on optional telemetry data files and the SPICE/Arrow/CSV ecosystem imported at load time. A missing manifest, unavailable GRAM SPICE directory, or solver failure is surfaced through the runner’s result/report path rather than repaired by the analysis layer. The strict defaults also assume reference samples use the same units, frame conventions, and timestep alignment as the simulation output.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification.jl` and its included `telemetry_verification/` implementation files.
