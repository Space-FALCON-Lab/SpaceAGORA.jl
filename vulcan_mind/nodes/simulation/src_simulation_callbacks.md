---
id: grp.src_simulation_callbacks
label: simulation/callbacks/
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
expands: module.simulation
tags:
- cluster
charts:
- simulation
origin: agent
---

# simulation/callbacks/

## Purpose
The integrator callbacks that surround the right-hand side: density and planet-frame staging, drag-state and event detection, guidance, navigation, control and thermal updates, data saving, and the caches that make GRAM affordable.

## Design & Implementation
`density_callbacks/` stages atmosphere per accepted step with batching, isolated GRAM pools and a vacuum-predicted spline cache; `gram_track_cache/` predicts density along a passage; `event_callbacks.jl` detects entry, exit, apsides and impact; `control_callbacks.jl` schedules thruster firings; `save_fields.jl` defines output columns; `registry.jl` holds runtime statistics.

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

- [[grp.src_simulation_engine|simulation/engine/]] · `members_out` → `members_in` · call · `src/simulation/engine/execution.jl:206-206`
- [[module.simulation|RuntimeServices]] · `api` → `members_in` · call · `src/simulation/callbacks/density_callbacks/assembly.jl`

**Downstream**

- `members_out` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:151-151`
- `members_out` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:112-112`
- `members_out` → [[core.reference_system_latlongtoned|latlongtoNED]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:37-37`
- `members_out` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:93-93`
- `members_out` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/save_fields.jl:94-94`
- `members_out` → [[core.runtime_types_callbackenvconfig|CallbackEnvConfig]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:178-178`
- `members_out` → [[core.runtime_types_gramtrackcacheconfig|GramTrackCacheConfig]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:145-145`
- `members_out` → [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:321-321`
- `members_out` → [[environment.get_density|getDensity]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:201-201`
- `members_out` → [[environment.simple_ephemerides_ephemerides_time_seconds|ephemerides_time_seconds]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:86-86`
- `members_out` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/callbacks/registry.jl:78-78`
- `members_out` → [[gnc.momentum_manager_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `members_out` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:148-148`
- `members_out` → [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `members_out` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/callbacks/registry.jl:78-78`
- `members_out` → [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `members_out` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/callbacks/registry.jl:78-78`
- `members_out` → [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `members_out` → [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `members_out` → [[gncz.navigation_hooks_calcnavigationeffect_bang|calcNavigationEffect!]] · `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:10-10`
- `members_out` → [[grp.cli|SpaceAGORACLI — internals]] · `members_in` · call · `src/simulation/callbacks/event_callbacks.jl:148-148`
- `members_out` → [[grp.src_analysis_verification|analysis/verification/]] · `members_in` · call · `src/simulation/callbacks/event_callbacks.jl:151-151`
- `members_out` → [[grp.src_core_interfaces|core/interfaces/]] · `members_in` · call · `src/simulation/callbacks/save_fields.jl:93-93`
- `members_out` → [[grp.src_environment_atmosphere|environment/atmosphere/]] · `members_in` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:321-321`
- `members_out` → [[grp.src_environment_ephemerides|environment/ephemerides/]] · `members_in` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:86-86`
- `members_out` → [[grp.src_gnc_control|gnc/control/]] · `members_in` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- `members_out` → [[grp.src_gnc_guidance|gnc/guidance/]] · `members_in` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:148-148`
- `members_out` → [[grp.src_parallel_policy|parallel/policy/]] · `members_in` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:36-36`
- `members_out` → [[grp.src_simulation_engine|simulation/engine/]] · `members_in` · call · `src/simulation/callbacks/density_callbacks/config.jl:303-303`
<!-- vulcan:connections:end -->

## Limitations
Callbacks fire on accepted steps only, so intermediate RK stages read whatever the last callback staged.

## Provenance
Macro block generated by `vulcan compile` (D13) from `src/simulation/callbacks`; groups 183 nodes.
Its members are listed in chart `simulation-simulation-callbacks`.
