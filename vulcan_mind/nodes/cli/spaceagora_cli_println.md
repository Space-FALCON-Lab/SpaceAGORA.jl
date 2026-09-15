---
id: cli.spaceagora_cli_println
label: println
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: println
  lines:
  - 61
  - 61
inputs:
- id: io
  type: Any
  units: n/a
  required: true
  description: Positional argument `io`.
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
  description: Return value of `println`. Returns `$(v)")`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
charts:
- cli
origin: agent
---

# println

## Purpose
This node was extracted from the `println(io, "cmd=$(full)")` call on line 61 inside `_run_subprocess`. It marks the final line of the `--print-only` dry-run report, where the fully interpolated Julia command is echoed instead of being executed.

## Design & Implementation
Within `_run_subprocess`, when `print_only` is true the function prints `project=<DOT_AGORA_PROJECT>`, `script=<script>`, an optional `env:` block listing each `k=v` pair from `env_pairs`, and lastly `cmd=<full>` where `full` is the `Cmd` object `` `$cmd --project=$DOT_AGORA_PROJECT $script $script_args` ``. Printing a `Cmd` uses Julia's backtick display, so arguments containing spaces appear quoted. The function then returns 0 without running anything or creating the output directory.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `io` | Any | n/a | yes | Positional argument `io`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `println`. Returns `$(v)")`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support_run_and_report|run_and_report]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:184-184`
- [[analysis.manifest_parsing__parse_maneuver_config|_parse_maneuver_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:264-264`
- [[analysis.runner__final_run_or_reused_eval|_final_run_or_reused_eval]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:328-328`
- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:361-361`
- [[analysis.scenario_builders__make_tabulated_flight_density_model|_make_tabulated_flight_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:239-239`
- [[analysis.scenario_builders__make_time_tabulated_density_model|_make_time_tabulated_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:262-262`
- [[analysis.scenario_builders__scenario_initial_condition|_scenario_initial_condition]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:531-531`
- [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:423-423`
- [[analysis.telemetry_loading__load_time_aligned_telemetry|_load_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:345-345`
- [[cli.assets_render_asset_manifest|render_asset_manifest]] · `callees` → `callers` · call · `src/cli/assets.jl:120-120`
- [[cli.assets_render_asset_report|render_asset_report]] · `callees` → `callers` · call · `src/cli/assets.jl:107-107`
- [[cli.assets_setup_open_assets|setup_open_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:138-138`
- [[cli.spaceagora_cli__print_usage|_print_usage]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:41-41`
- [[cli.spaceagora_cli__run_subprocess|_run_subprocess]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:56-56`
- [[core.reference_system_latlongtooe|latlongtoOE]] · `callees` → `callers` · call · `src/core/interfaces/reference_system.jl:306-306`
- [[dynamics.aerodynamic_wrench_models_compute_link_wrench_bang|compute_link_wrench!]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:911-911`
- [[envana.ana_calibration_estimate_event_biases|_estimate_event_biases]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:47-47`
- [[envana.ana_ic_fit_fit_initial_state|fit_initial_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:178-178`
- [[gnc.eom_predictor_shooting_residual_bang|shooting_residual!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1041-1041`
- [[gnc.planner_comparison__rpo_comparison_progress_line_bang|_rpo_comparison_progress_line!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:88-88`
- [[gnc.planner_comparison_rpo_write_namedtuple_csv|rpo_write_namedtuple_csv]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:663-663`
- [[gnc.propulsive_maneuvers__trace_maneuver_event_bang|_trace_maneuver_event!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:303-303`
- [[gnc.targeting_solver__target_planning_impl|_target_planning_impl]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:24-24`
- [[gnc.targeting_solver_func_targeting_heatload|func_targeting_heatload]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:154-154`
- [[gnc.targeting_solver_func_targeting_num_int|func_targeting_num_int]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:128-128`
- [[gnc.tracking_executor_df|df]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:155-155`
- [[gnc.tracking_executor_f|f]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:60-60`
- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:204-204`
- [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:155-155`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1041-1041`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:204-204`
- [[grp.src_analysis_verification|analysis/verification/]] · `members_out` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:184-184`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:303-303`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1041-1041`
- [[grp.src_parallel_routing|parallel/routing/]] · `members_out` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:603-603`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:112-112`
- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `callers` · call · `src/simulation/engine/reporting.jl:31-31`
- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:603-603`
- [[simulation.dynamics_rhs__accumulate_control_effectors_bang|_accumulate_control_effectors!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:199-199`
- [[simulation.event_callbacks_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:56-56`
- [[simulation.event_callbacks_affect_downcrossing_bang|affect_downcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:16-16`
- [[simulation.event_callbacks_affect_upcrossing_bang|affect_upcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:148-148`
- [[simulation.event_callbacks_condition|condition]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:223-223`
- [[simulation.event_callbacks_get_entry_end_callback|get_entry_end_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:112-112`
- [[simulation.reporting__debug_print_nan_parameter_paths_bang|_debug_print_nan_parameter_paths!]] · `callees` → `callers` · call · `src/simulation/engine/reporting.jl:31-31`
- [[simulation.rhs_calibration__run_rhs_sweep_bang|_run_rhs_sweep!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:261-261`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1848-1848`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1899-1899`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1794-1794`
- [[simulation_a.event_callbacks_get_impact_callback|get_impact_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:16-16`
- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:179-179`
- [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callees` → `callers` · call · `src/simulation/engine/rhs_calibration.jl:333-333`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The echoed command does not include the environment overrides inline, so copying `cmd=` into a shell does not reproduce a run that depends on `SPACEAGORA_*` variables; the user must also export the `env:` lines. `Base.julia_cmd()` embeds the current process's flags, which may differ from what a user would type.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 61.
