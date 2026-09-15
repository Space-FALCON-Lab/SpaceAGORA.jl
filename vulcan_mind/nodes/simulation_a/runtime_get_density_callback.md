---
id: simulation_a.runtime_get_density_callback
label: get_density_callback
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/runtime.jl
  symbol: get_density_callback
  lines:
  - 212
  - 356
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: num_sats_effectors_args
  type: Tuple{Int, Tuple, SimulationConfiguration}
  units: n/a
  required: true
  description: Spacecraft count, the dynamic effector tuple used to detect a J2 gravity
    effector, and the run configuration.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: density_callback
  type: DiscreteCallback
  units: n/a
  description: Per-step callback that writes density, temperature and wind for every
    active spacecraft into the shared buffers consumed by the aerodynamic and thermal
    paths.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# get_density_callback

## Purpose
`get_density_callback` builds the callback that samples the atmosphere once per accepted solver step for every spacecraft and publishes the result into shared buffers. Aerodynamic effectors and the thermal callback read those buffers rather than querying the atmosphere themselves, so this callback is the single point where the cost of an atmosphere model is paid.

## Model & Assumptions
Four evaluation strategies are selected at run time. The per-spacecraft path calls `_density_state_from_kinematics!`, which may consult the GRAM track cache before falling back to a direct query. The batch path packs altitude, latitude and longitude into preallocated shared arrays and issues one `getDensityBatch!` call. The isolated-pool path replaces that call with `_gram_isolated_pool_batch_eval!` when a threaded pool of private GRAM copies is available. Batching is only chosen when the active model can evaluate all spacecraft together without per-satellite cache state, which is why `_gram_track_cache_enabled` disqualifies it. J2-aware track-cache targeting is enabled only when a J2 gravity effector is actually present in the effector tuple.

## Design & Implementation
Two methods are exposed: a two-argument convenience form that pulls the effector tuple out of `args.dynamics_model`, and the three-argument form that does the work. GRAM knobs are read from the run-scoped snapshot through `_callback_env_config(p)` at invocation time rather than captured from live ENV at construction, so the callback and the RHS-side sampler always agree. The kinematics gather loop is itself threaded when `ParallelPolicy` approves, and elapsed nanoseconds are reported back through `record_policy_observation!(:density_callback, ...)`. The callback's initializer runs the affect immediately so buffers are populated before the first derivative call.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `num_sats_effectors_args` | Tuple{Int, Tuple, SimulationConfiguration} | n/a | yes | Spacecraft count, the dynamic effector tuple used to detect a J2 gravity effector, and the run configuration. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `density_callback` | DiscreteCallback | n/a | — | Per-step callback that writes density, temperature and wind for every active spacecraft into the shared buffers consumed by the aerodynamic and thermal paths. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:161-161`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`
- `callees` → [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:321-321`
- `callees` → [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:291-291`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:345-345`
- `callees` → [[simulation.assembly__uses_j2_gravity_effector|_uses_j2_gravity_effector]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:221-221`
- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:223-223`
- `callees` → [[simulation.config__density_batch_enabled|_density_batch_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:260-260`
- `callees` → [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:255-255`
- `callees` → [[simulation.event_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:249-249`
- `callees` → [[simulation.event_callbacks_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:247-247`
- `callees` → [[simulation.interpolation__gram_track_cache_enabled|_gram_track_cache_enabled]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:262-262`
- `callees` → [[simulation.model_selection__density_batch_model_for_callback|_density_batch_model_for_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:261-261`
- `callees` → [[simulation.model_selection__density_model_for_sat|_density_model_for_sat]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:224-224`
- `callees` → [[simulation.model_selection__gram_isolated_pool_batch_model_for_callback|_gram_isolated_pool_batch_model_for_callback]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:267-267`
- `callees` → [[simulation.planet_frame_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:249-249`
- `callees` → [[simulation.planet_frame_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:247-247`
- `callees` → [[simulation.registry__gram_runtime_stats_update_bang|_gram_runtime_stats_update!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:282-282`
- `callees` → [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:227-227`
- `callees` → [[simulation.runtime__extract_mass_kg|_extract_mass_kg]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:232-232`
- `callees` → [[simulation.runtime__stage_environment_kinematics|_stage_environment_kinematics]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:226-226`
- `callees` → [[simulation.runtime__write_density_buffers_bang|_write_density_buffers!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:243-243`
- `callees` → [[simulation.runtime__write_density_time_buffers_bang|_write_density_time_buffers!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:334-334`
- `callees` → [[simulation.runtime_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:249-249`
- `callees` → [[simulation.runtime_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:247-247`
- `callees` → [[simulation.runtime_update_density_sat_bang|update_density_sat!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:222-222`
- `callees` → [[simulation.thermal_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:249-249`
- `callees` → [[simulation.thermal_callbacks_condition|condition]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:247-247`
- `callees` → [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:307-307`
<!-- vulcan:connections:end -->

## Limitations
Threading correctness rests on distinct spacecraft indices touching disjoint buffer slots and on the density model being safe to call concurrently, which is why the isolated pool exists for GRAM. The condition is unconditionally true, so cost scales with accepted step count and spacecraft count. Buffers are overwritten in place, so any consumer needing a history must copy the values within the same step.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/runtime.jl:211-356`.
