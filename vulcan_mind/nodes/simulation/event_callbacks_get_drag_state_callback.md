---
id: simulation.event_callbacks_get_drag_state_callback
label: get_drag_state_callback
kind: function
source:
  file: src/simulation/callbacks/event_callbacks.jl
  symbol: get_drag_state_callback
  lines:
  - 138
  - 138
inputs:
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  type: VectorContinuousCallback
  units: n/a
  description: Return value of `get_drag_state_callback`. Returns `VectorContinuousCallback(condition!,
    affect_upcrossing!, affect_downcrossing!, n`.
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

# get_drag_state_callback

## Purpose
Builds the `VectorContinuousCallback` that tracks whether each satellite is inside the sensible atmosphere and switches integrator settings accordingly. Crossing the entry interface downward tightens step size and tolerances for drag-pass integration; crossing upward relaxes them for orbital coasting and invalidates the vacuum GRAM cache.

## Design & Implementation
Takes `num_sats::Int` and returns `VectorContinuousCallback(condition!, affect_upcrossing!, affect_downcrossing!, num_sats)`. `condition!` computes `alt - EI*1e3` per satellite, positive above the atmosphere. Both affects write `p.shared_buffers.in_atmosphere[idx]` and stamp `in_atmosphere_sample_t[idx] = Float64(integrator.t)`. The downcrossing sets `integrator.opts.dtmax = dt_max_atmosphere` and applies `_callback_tolerances_for_phase(reltol, abstol, args, true)`; the upcrossing sets `dt_max_orbit`, applies the `false` phase tolerances, marks `vacuum_gram_caches[idx].valid = false` when present, and calls `schedule_event_driven_thruster_controls!(integrator, idx)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VectorContinuousCallback | n/a | — | Return value of `get_drag_state_callback`. Returns `VectorContinuousCallback(condition!, affect_upcrossing!, affect_downcrossing!, n`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.assembly_get_callbacks|get_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/assembly.jl:177-177`

**Downstream**

- `callees` → [[simulation.event_callbacks_condition_bang|condition!]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:139-139`
- `callees` → [[simulation.registry__simulation_engine_module|_simulation_engine_module]] · `callers` · call · `src/simulation/callbacks/event_callbacks.jl:141-141`
<!-- vulcan:connections:end -->

## Limitations
The callback mutates `integrator.opts` globally for all satellites, so in a multi-satellite run the last crossing wins: one satellite entering the atmosphere forces atmospheric step limits on every satellite in the shared state. `in_atmosphere` is never initialised here; the initial flag must be set correctly by the caller or the first crossing is misinterpreted. Cache invalidation is skipped silently when `idx > length(vacuum_gram_caches)`. Altitude uses the equatorial radius only.

## Provenance
Mapped from `src/simulation/callbacks/event_callbacks.jl` line 138.
