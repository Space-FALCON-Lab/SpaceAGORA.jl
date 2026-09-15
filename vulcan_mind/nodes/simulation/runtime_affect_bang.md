---
id: simulation.runtime_affect_bang
label: affect!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: affect!
  lines:
  - 249
  - 249
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `affect!`; mutates `integrator` in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# affect!

## Purpose
The body of the density callback: decide whether to evaluate all satellites in one batch, in parallel, or serially, run the sampling, and record timing for the adaptive thread policy.

## Design & Implementation
Reads density models and the run-scoped env config, then asks `_density_callback_thread_decision` for a threading verdict. It prefers a batch model when `_density_batch_enabled` finds one shared model without track-cache state, falling back to the GRAM isolated-pool batch model; a GRAM batch model also sets `use_gram_isolated_pool`. In batch mode it fills the altitude, latitude and longitude scratch buffers — threaded or not — then calls the pooled evaluator or `getDensityBatch!`, and stamps sample times with `_write_density_time_buffers!`. Otherwise it runs `update_density_sat!` per satellite, threaded via `threaded_foreach_persistent` with the decided allotment or in a plain loop. If the policy applied, elapsed nanoseconds are reported through `record_policy_observation!`. The callback is also invoked once from `initialize` so buffers are valid before the first RHS call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `affect!`; mutates `integrator` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.event_callbacks_condition|condition]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:205-205`
- [[simulation.planet_frame_init_affect_bang|init_affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:87-87`
- [[simulation_a.planet_frame_update_planet_frame_callback|update_planet_frame_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:71-71`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:249-249`
- [[simulation_a.thermal_callbacks_get_thermal_callback|get_thermal_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:68-68`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:293-293`
- `callees` → [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:321-321`
- `callees` → [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:291-291`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:345-345`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:254-254`
- `callees` → [[simulation.config__density_batch_enabled|_density_batch_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:260-260`
- `callees` → [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:255-255`
- `callees` → [[simulation.interpolation__gram_track_cache_enabled|_gram_track_cache_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:262-262`
- `callees` → [[simulation.model_selection__density_batch_model_for_callback|_density_batch_model_for_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:261-261`
- `callees` → [[simulation.model_selection__gram_isolated_pool_batch_model_for_callback|_gram_isolated_pool_batch_model_for_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:267-267`
- `callees` → [[simulation.registry__gram_runtime_stats_update_bang|_gram_runtime_stats_update!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:282-282`
- `callees` → [[simulation.runtime__stage_environment_kinematics|_stage_environment_kinematics]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:293-293`
- `callees` → [[simulation.runtime__write_density_time_buffers_bang|_write_density_time_buffers!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:334-334`
- `callees` → [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:337-337`
- `callees` → [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:307-307`
<!-- vulcan:connections:end -->

## Limitations
Batch mode is chosen only when no satellite uses the track cache, so a mixed constellation silently degrades to per-satellite evaluation; the thread decision is recomputed on every step, and its cost is not excluded from the elapsed time it reports back to the policy.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl` line 249.
