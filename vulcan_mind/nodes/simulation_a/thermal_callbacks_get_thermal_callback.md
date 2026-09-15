---
id: simulation_a.thermal_callbacks_get_thermal_callback
label: get_thermal_callback
kind: function
source:
  file: src/simulation/callbacks/thermal_callbacks.jl
  symbol: get_thermal_callback
  lines:
  - 60
  - 94
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: num_sats_and_config
  type: Tuple{Int, SimulationConfiguration}
  units: n/a
  required: true
  description: Spacecraft count and the run configuration, used to size the per-spacecraft
    loop and to drive the parallel-policy decision.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: thermal_callback
  type: DiscreteCallback
  units: n/a
  description: Per-step callback that writes stage heat rates for every spacecraft
    into the shared thermal buffers integrated as accumulated heat load.
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
# get_thermal_callback

## Purpose
`get_thermal_callback` builds the callback that evaluates aerothermal heating for every spacecraft once per accepted solver step. It is installed only when `_requires_thermal_callback` finds a thermal model attached to the effector set, so vacuum and non-heating runs never pay for it.

## Model & Assumptions
Heating is computed per stage through `_compute_stage_heat_rates!`, which calls `getHeatRate` on the vehicle thermal model. The callback runs with `use_buffered_density=true`, meaning it consumes the density, temperature and wind values the density callback has already written into shared buffers for this step rather than querying the atmosphere itself. That coupling is why `get_callbacks` installs the density callback before the thermal one. Computed rates are floored at zero and non-finite values are replaced by zero, so a model returning `NaN` degrades to no heating rather than poisoning the integrated heat load.

## Design & Implementation
The structure mirrors the density and control callbacks: an unconditionally true condition, an affect that asks `_thermal_callback_thread_decision` whether to thread, and a `ParallelPolicy.threaded_foreach_persistent(:thermal_callback, ...)` loop when it approves, falling back to a plain `@inbounds` loop otherwise. Elapsed nanoseconds are reported through `record_policy_observation!` whenever the policy was actually applied, feeding the adaptive threshold. The callback's initializer invokes the affect directly so heat rates are populated before the first derivative evaluation instead of remaining zero through the opening step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `num_sats_and_config` | Tuple{Int, SimulationConfiguration} | n/a | yes | Spacecraft count and the run configuration, used to size the per-spacecraft loop and to drive the parallel-policy decision. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `thermal_callback` | DiscreteCallback | n/a | — | Per-step callback that writes stage heat rates for every spacecraft into the shared thermal buffers integrated as accumulated heat load. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:165-165`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:76-76`
- `callees` → [[parallel.thread_execution_threaded_foreach_persistent|threaded_foreach_persistent]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:75-75`
- `callees` → [[parcore.observation_tracking_record_policy_observation_bang|record_policy_observation!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:84-84`
- `callees` → [[simulation.config__thermal_callback_thread_decision|_thermal_callback_thread_decision]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:71-71`
- `callees` → [[simulation.event_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:68-68`
- `callees` → [[simulation.event_callbacks_condition|condition]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:66-66`
- `callees` → [[simulation.planet_frame_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:68-68`
- `callees` → [[simulation.planet_frame_condition|condition]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:66-66`
- `callees` → [[simulation.runtime_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:68-68`
- `callees` → [[simulation.runtime_condition|condition]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:66-66`
- `callees` → [[simulation.thermal_callbacks__compute_stage_heat_rates_bang|_compute_stage_heat_rates!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:62-62`
- `callees` → [[simulation.thermal_callbacks_affect_bang|affect!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:68-68`
- `callees` → [[simulation.thermal_callbacks_condition|condition]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:66-66`
- `callees` → [[simulation.thermal_callbacks_update_thermal_sat_bang|update_thermal_sat!]] · `callers` · call · `src/simulation/callbacks/thermal_callbacks.jl:61-61`
<!-- vulcan:connections:end -->

## Limitations
Threaded evaluation requires the thermal model to be safe under concurrent calls on distinct spacecraft indices; models holding shared scratch arrays must be forced onto the serial path through the parallel-mode setting. Because the callback reads buffered density, any step in which the density callback did not run leaves it working from stale atmosphere values. The callback is omitted in gravity-backbone split mode.

## Provenance
Mapped from `src/simulation/callbacks/thermal_callbacks.jl:59-94`.
