---
id: simulation.event_callbacks_get_entry_end_callback
label: get_entry_end_callback
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: get_entry_end_callback
  lines:
  - 84
  - 84
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
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
  type: Union{Nothing, VectorContinuousCallback}
  units: n/a
  description: Return value of `get_entry_end_callback`. Returns `nothing` or `VectorContinuousCallback(condition!,
    nothing, affect_downcrossing!, num_sats)`.
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

# get_entry_end_callback

## Purpose
Builds a `VectorContinuousCallback` that counts atmospheric entry-interface crossings per satellite and terminates the simulation once every active satellite has recorded `SPACEAGORA_ENTRY_TARGET_COUNT` entries. It is intended for multi-pass aerobraking studies where the stop criterion is a number of drag passes rather than orbits or time.

## Design & Implementation
Arguments are `num_sats::Int` and `args::SimulationConfiguration`. `_entry_target_count()` reads the environment variable and the function throws `ArgumentError` when the target is not positive. The interface altitude is `args.environment_model.EI * 1e3` (EI in km converted to m). A closure-local `entry_counter = zeros(Int64, num_sats)` is captured. `condition!` writes `alt - entry_interface_m` for active satellites and `1.0` for inactive ones so they never trigger. `affect_downcrossing!` increments `entry_counter[idx]`, and when `completed_entries >= target_entries` scans all active satellites; if all have met the target it calls `terminate!` guarded by `applicable`. Returns `VectorContinuousCallback(condition!, nothing, affect_downcrossing!, num_sats)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, VectorContinuousCallback} | n/a | — | Return value of `get_entry_end_callback`. Returns `nothing` or `VectorContinuousCallback(condition!, nothing, affect_downcrossing!, num_sats)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:173-173`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:112-112`
- `callees` → [[simulation.assembly__entry_target_count|_entry_target_count]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:85-85`
- `callees` → [[simulation.event_callbacks_affect_downcrossing_bang|affect_downcrossing!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:103-103`
- `callees` → [[simulation.event_callbacks_condition_bang|condition!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:90-90`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:98-98`
- `callees` → [[simulation.registry_callback_verbose|callback_verbose]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:111-111`
<!-- vulcan:connections:end -->

## Limitations
The entry counter lives in the closure, not in `integrator.p`, so it is invisible to telemetry and is not reset if the same callback object is reused across `solve` calls. Only downcrossings count, so a satellite starting below the interface records its first entry on the second pass. The target is read from an environment variable at construction time, not from `args`, which makes the configuration implicit and non-reproducible from the configuration file alone. Altitude is spherical.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 84.
