---
id: simx.engine_execution_run_simulation
label: run_simulation
kind: function
source:
  file: src/simulation/engine/execution.jl
  symbol: run_simulation
  lines:
  - 143
  - 468
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: engine_api
  type: Module
  units: n/a
  required: true
  description: RuntimeServices namespace and the SimulationEngine/SimulationCampaigns
    scope it aggregates, supplying the SPICE lock and the SimulationModel types used
    here.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: 'Complete typed configuration: dynamics model, environment model, mission
    configuration, simulation settings and optional solver config.'
- id: isolate_state
  type: Bool
  units: n/a
  required: true
  description: When true the configuration is deepcopied so repeated or concurrent
    runs do not alias shared mutable model state.
- id: solver_cache
  type: Union{Nothing,SolverIntegratorCache}
  units: n/a
  required: true
  description: Optional reusable integrator cache that lets repeated solves of the
    same shape skip integrator allocation.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: solution
  type: Union{Nothing,ODESolution,NamedTuple}
  units: n/a
  description: Nothing by default; the final ODESolution when return_solution is set;
    a named tuple adding solver_mode, solver_trace, parallel_policy and spice_counters
    when return_solver_metadata is also set.
- id: result_artifacts
  type: String
  units: n/a
  description: Paths of the results CSV and optional results bundle written by the
    persistence layer.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simx
origin: agent
---
# run_simulation

## Purpose
`run_simulation` is the engine's single entry point for propagating one configuration. It validates the request, builds initial conditions and ODE parameters, allocates every workspace buffer, selects and runs the solver, and hands the saved history to the persistence layer.

## Model & Assumptions
The typed pipeline is SI-native in metres, seconds and kilograms. `_enforce_typed_normalize_policy!` therefore rejects `simulation_settings.normalize=true` outright unless the transition escape hatch is set, in which case it warns once. Three further guards run before any allocation: `_validate_orientation_inertia!`, `_validate_thermal_model_support!` and `_validate_ephemerides_support!`, plus `_warn_density_without_atmospheric_effector` for configurations that request density with nothing to consume it.

## Design & Implementation
The whole body executes inside `SimulationModel.ParallelPolicy.with_policy_context()`, so thread-budget decisions are scoped to the run. The effective solver configuration is `args.solver_config` when present and `_solver_config_from_env()` otherwise. After `build_initial_conditions`, an `ODEParams` is constructed for `n_sats` spacecraft and then twelve initializers run in a fixed order: in-atmosphere flags, the runtime env-config snapshot, heat-rate buffers, save-cache buffers, density model instances and cache buffers, the GRAM isolated pool, harmonics, n-body and aero workspaces, and the n-body, SRP-sun and planet-frame ephemeris cache buffers. Integration is either a checkpointed segment loop that writes a checkpoint after each segment via `_write_checkpoint!`, or a single span solved by `_solve_with_solver_policy`. Every failure path calls `_try_save_simulation_results_if_enabled!` before rethrowing, so partial history survives a crash.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `engine_api` | Module | n/a | yes | RuntimeServices namespace and the SimulationEngine/SimulationCampaigns scope it aggregates, supplying the SPICE lock and the SimulationModel types used here. |
| in | `args` | SimulationConfiguration | n/a | yes | Complete typed configuration: dynamics model, environment model, mission configuration, simulation settings and optional solver config. |
| in | `isolate_state` | Bool | n/a | yes | When true the configuration is deepcopied so repeated or concurrent runs do not alias shared mutable model state. |
| in | `solver_cache` | Union{Nothing,SolverIntegratorCache} | n/a | yes | Optional reusable integrator cache that lets repeated solves of the same shape skip integrator allocation. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `solution` | Union{Nothing,ODESolution,NamedTuple} | n/a | — | Nothing by default; the final ODESolution when return_solution is set; a named tuple adding solver_mode, solver_trace, parallel_policy and spice_counters when return_solver_metadata is also set. |
| out | `result_artifacts` | String | n/a | — | Paths of the results CSV and optional results bundle written by the persistence layer. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support_run_and_report|run_and_report]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:179-179`
- [[analysis.runner__run_once|_run_once]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:47-47`
- [[misc.precompile_workload_run_spaceagora_precompile_workload|_run_spaceagora_precompile_workload]] · `callees` → `callers` · call · `src/precompile_workload.jl:47-47`
- [[simulation.monte_carlo__run_monte_carlo_process|_run_monte_carlo_process]] · `callees` → `callers` · call · `src/simulation/campaigns/monte_carlo.jl:223-223`
- [[simx.campaigns_constellation_ensemble_run_constellation_ensemble|run_constellation_ensemble]] · `callees` → `callers` · call · `src/simulation/campaigns/constellation_ensemble.jl:160-160`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:492-492`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/engine/execution.jl:385-385`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/execution.jl:179-179`
- `callees` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/simulation/engine/execution.jl:219-219`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/simulation/engine/execution.jl:221-221`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/execution.jl:357-357`
- `callees` → [[io.io_serialization__load_checkpoint|_load_checkpoint]] · `callers` · call · `src/simulation/engine/execution.jl:237-237`
- `callees` → [[misc.io_serialization_write_checkpoint__write_checkpoint_bang|_write_checkpoint!]] · `callers` · call · `src/simulation/engine/execution.jl:387-387`
- `callees` → [[parallel.env_config_persistent_hints_state_reset_requested|persistent_hints_state_reset_requested]] · `callers` · call · `src/simulation/engine/execution.jl:170-170`
- `callees` → [[parallel.persistent_hints_reset_persistent_hint_state_bang|reset_persistent_hint_state!]] · `callers` · call · `src/simulation/engine/execution.jl:171-171`
- `callees` → [[parallel.policy_telemetry_reset_policy_telemetry_bang|reset_policy_telemetry!]] · `callers` · call · `src/simulation/engine/execution.jl:169-169`
- `callees` → [[parcore.context_with_policy_context|with_policy_context]] · `callers` · call · `src/simulation/engine/execution.jl:152-152`
- `callees` → [[parcore.policy_telemetry_policy_telemetry_snapshot|policy_telemetry_snapshot]] · `callers` · call · `src/simulation/engine/execution.jl:452-452`
- `callees` → [[parcore.runtime_types_odeparams|ODEParams]] · `callers` · call · `src/simulation/engine/execution.jl:184-184`
- `callees` → [[simulation.dynamics_rhs_build_initial_conditions|build_initial_conditions]] · `callers` · call · `src/simulation/engine/execution.jl:177-177`
- `callees` → [[simulation.execution__append_backbone_saved_segment_bang|_append_backbone_saved_segment!]] · `callers` · call · `src/simulation/engine/execution.jl:359-359`
- `callees` → [[simulation.execution__append_checkpoint_saved_segment_bang|_append_checkpoint_saved_segment!]] · `callers` · call · `src/simulation/engine/execution.jl:361-361`
- `callees` → [[simulation.execution__build_block_diagonal_jac_prototype|_build_block_diagonal_jac_prototype]] · `callers` · call · `src/simulation/engine/execution.jl:308-308`
- `callees` → [[simulation.execution__build_typed_solver_problem|_build_typed_solver_problem]] · `callers` · call · `src/simulation/engine/execution.jl:340-340`
- `callees` → [[simulation.execution__save_simulation_results_if_enabled_bang|_save_simulation_results_if_enabled!]] · `callers` · call · `src/simulation/engine/execution.jl:434-434`
- `callees` → [[simulation.execution__try_save_simulation_results_if_enabled_bang|_try_save_simulation_results_if_enabled!]] · `callers` · call · `src/simulation/engine/execution.jl:344-344`
- `callees` → [[simulation.from_env__engine_env_get|_engine_env_get]] · `callers` · call · `src/simulation/engine/execution.jl:204-204`
- `callees` → [[simulation.from_env__solver_config_from_env|_solver_config_from_env]] · `callers` · call · `src/simulation/engine/execution.jl:158-158`
- `callees` → [[simulation.reporting__debug_print_nan_parameter_paths_bang|_debug_print_nan_parameter_paths!]] · `callers` · call · `src/simulation/engine/execution.jl:281-281`
- `callees` → [[simulation.resume_checkpoint__load_checkpoint|_load_checkpoint]] · `callers` · call · `src/simulation/engine/execution.jl:237-237`
- `callees` → [[simulation.setup__initialize_aero_workspace_buffers_bang|_initialize_aero_workspace_buffers!]] · `callers` · call · `src/simulation/engine/execution.jl:197-197`
- `callees` → [[simulation.setup__initialize_density_cache_buffers_bang|_initialize_density_cache_buffers!]] · `callers` · call · `src/simulation/engine/execution.jl:193-193`
- `callees` → [[simulation.setup__initialize_density_model_instances_bang|_initialize_density_model_instances!]] · `callers` · call · `src/simulation/engine/execution.jl:192-192`
- `callees` → [[simulation.setup__initialize_gram_isolated_pool_buffers_bang|_initialize_gram_isolated_pool_buffers!]] · `callers` · call · `src/simulation/engine/execution.jl:194-194`
- `callees` → [[simulation.setup__initialize_harmonics_workspace_buffers_bang|_initialize_harmonics_workspace_buffers!]] · `callers` · call · `src/simulation/engine/execution.jl:195-195`
- `callees` → [[simulation.setup__initialize_heat_rate_buffers_bang|_initialize_heat_rate_buffers!]] · `callers` · call · `src/simulation/engine/execution.jl:190-190`
- `callees` → [[simulation.setup__initialize_in_atmosphere_flags_bang|_initialize_in_atmosphere_flags!]] · `callers` · call · `src/simulation/engine/execution.jl:185-185`
- `callees` → [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callers` · call · `src/simulation/engine/execution.jl:226-226`
- `callees` → [[simulation.setup__initialize_nbody_ephemeris_cache_buffer_bang|_initialize_nbody_ephemeris_cache_buffer!]] · `callers` · call · `src/simulation/engine/execution.jl:198-198`
- `callees` → [[simulation.setup__initialize_nbody_workspace_buffers_bang|_initialize_nbody_workspace_buffers!]] · `callers` · call · `src/simulation/engine/execution.jl:196-196`
- `callees` → [[simulation.setup__initialize_planet_frame_cache_buffer_bang|_initialize_planet_frame_cache_buffer!]] · `callers` · call · `src/simulation/engine/execution.jl:200-200`
- `callees` → [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callers` · call · `src/simulation/engine/execution.jl:228-228`
- `callees` → [[simulation.setup__initialize_runtime_env_config_bang|_initialize_runtime_env_config!]] · `callers` · call · `src/simulation/engine/execution.jl:189-189`
- `callees` → [[simulation.setup__initialize_save_cache_buffers_bang|_initialize_save_cache_buffers!]] · `callers` · call · `src/simulation/engine/execution.jl:191-191`
- `callees` → [[simulation.setup__initialize_spice_rhs_memo_mode_bang|_initialize_spice_rhs_memo_mode!]] · `callers` · call · `src/simulation/engine/execution.jl:201-201`
- `callees` → [[simulation.setup__initialize_srp_sun_cache_buffer_bang|_initialize_srp_sun_cache_buffer!]] · `callers` · call · `src/simulation/engine/execution.jl:199-199`
- `callees` → [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callers` · call · `src/simulation/engine/execution.jl:227-227`
- `callees` → [[simulation.setup__reset_spice_rhs_memo_bang|_reset_spice_rhs_memo!]] · `callers` · call · `src/simulation/engine/execution.jl:203-203`
- `callees` → [[simulation.setup__reset_spice_runtime_counters_bang|_reset_spice_runtime_counters!]] · `callers` · call · `src/simulation/engine/execution.jl:202-202`
- `callees` → [[simulation.setup__spice_runtime_counters_snapshot|_spice_runtime_counters_snapshot]] · `callers` · call · `src/simulation/engine/execution.jl:461-461`
- `callees` → [[simulation.setup__typed_checkpoint_enabled|_typed_checkpoint_enabled]] · `callers` · call · `src/simulation/engine/execution.jl:229-229`
- `callees` → [[simulation.setup__validate_ephemerides_support_bang|_validate_ephemerides_support!]] · `callers` · call · `src/simulation/engine/execution.jl:166-166`
- `callees` → [[simulation.setup__validate_orientation_inertia_bang|_validate_orientation_inertia!]] · `callers` · call · `src/simulation/engine/execution.jl:164-164`
- `callees` → [[simulation.setup__validate_thermal_model_support_bang|_validate_thermal_model_support!]] · `callers` · call · `src/simulation/engine/execution.jl:165-165`
- `callees` → [[simulation.setup__warn_density_without_atmospheric_effector|_warn_density_without_atmospheric_effector]] · `callers` · call · `src/simulation/engine/execution.jl:167-167`
- `callees` → [[simulation.solver_policy__build_solver_tolerances|_build_solver_tolerances]] · `callers` · call · `src/simulation/engine/execution.jl:302-302`
- `callees` → [[simulation.solver_policy__gravity_backbone_time_reached|_gravity_backbone_time_reached]] · `callers` · call · `src/simulation/engine/execution.jl:388-388`
- `callees` → [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callers` · call · `src/simulation/engine/execution.jl:210-210`
- `callees` → [[simulation_a.save_fields_default_save_fields|default_save_fields]] · `callers` · call · `src/simulation/engine/execution.jl:206-206`
- `callees` → [[simx.engine_reporting_enforce_typed_normalize_policy__enforce_typed_normalize_policy_bang|_enforce_typed_normalize_policy!]] · `callers` · call · `src/simulation/engine/execution.jl:163-163`
- `callees` → [[simx.engine_resume_checkpoint_write_checkpoint__write_checkpoint_bang|_write_checkpoint!]] · `callers` · call · `src/simulation/engine/execution.jl:387-387`
- `callees` → [[simx.engine_rhs_calibration_calibrate_rhs_plan_if_needed__calibrate_rhs_plan_if_needed_bang|_calibrate_rhs_plan_if_needed!]] · `callers` · call · `src/simulation/engine/execution.jl:313-313`
- `callees` → [[simx.engine_solver_policy_solve_with_solver_policy|_solve_with_solver_policy]] · `callers` · call · `src/simulation/engine/execution.jl:342-342`
<!-- vulcan:connections:end -->

## Limitations
A non-successful `retcode` is converted into a plain `ErrorException` carrying the code, losing the solver's own diagnostic object. With checkpointing active and `return_solution=true` the caller receives only the final segment's `ODESolution` and is warned that the history is not stitched. The telemetry reset and the parallel-policy snapshot are wrapped in bare `try`/`catch` blocks that swallow any error, so a broken telemetry backend degrades silently.

## Provenance
Mapped from `src/simulation/engine/execution.jl:143-468`, with the block-diagonal Jacobian prototype builder at line 1 and the save-on-failure wrapper at line 134 of the same file.
