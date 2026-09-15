---
id: simulation.thermal_callbacks_affect_bang
label: affect!
kind: function
source:
  file: src/simulation/callbacks/thermal_callbacks.jl
  symbol: affect!
  lines:
  - 68
  - 68
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
Body of the thermal `DiscreteCallback`. It updates the heat rates of all `num_sats` spacecraft at the current integrator time, either serially or across threads, and records a timing observation for the adaptive parallel policy.

## Design & Implementation
It pulls `p = integrator.p` and `u = integrator.u`, then asks `_thermal_callback_thread_decision(p, num_sats)` whether to thread and with what allotment. The threaded branch uses `ParallelPolicy.threaded_foreach_persistent(:thermal_callback, num_sats, decision.allotment)`, which reuses a persistent worker set rather than spawning tasks per firing; the serial branch is a plain `@inbounds` loop. Both call `update_thermal_sat!(i, p, u, Float64(integrator.t))`. When `decision.policy_applied` is set, elapsed nanoseconds from `time_ns()` are fed back through `ParallelPolicy.record_policy_observation!` so future decisions adapt. The same function is wired as the callback's `initialize` hook, so heat rates are populated before the first step.

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

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:76-76`
- `callees` → [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:75-75`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:84-84`
- `callees` → [[simulation.config__thermal_callback_thread_decision|_thermal_callback_thread_decision]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:71-71`
- `callees` → [[simulation.thermal_callbacks_update_thermal_sat_bang|update_thermal_sat!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:76-76`
<!-- vulcan:connections:end -->

## Limitations
Thread safety rests on each spacecraft index owning a disjoint slice of `p.shared_buffers.heat_rates`; any effector that writes across spacecraft from inside the heat-rate path would race. Timing is wall-clock `time_ns()`, so measurements taken under external load bias the policy toward or away from threading. Errors raised in one worker abort the whole callback, and partial updates already written to the buffers are left in place.

## Provenance
Mapped from `src/simulation/callbacks/thermal_callbacks.jl` line 68.
