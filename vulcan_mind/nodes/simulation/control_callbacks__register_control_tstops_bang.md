---
id: simulation.control_callbacks__register_control_tstops_bang
label: _register_control_tstops!
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: _register_control_tstops!
  lines:
  - 16
  - 16
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
- id: control_model
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `control_model`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: Nothing
  units: n/a
  description: Return value of `_register_control_tstops!`; mutates `integrator` in
    place. Returns `nothing`.
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

# _register_control_tstops!

## Purpose
Registers both ends of a thruster burn for spacecraft `sat_idx` as integrator stop points, so the solver lands exactly on the ignition and cutoff instants instead of stepping across the thrust discontinuity.

## Design & Implementation
Marked `@inline`. It bounds-checks `sat_idx` against `length(control_model.start_burn_time)` and returns early when the index is outside `1:n`, which makes the routine safe to call for spacecraft that the model does not cover. It then passes `control_model.start_burn_time[sat_idx]` and `control_model.stop_burn_time[sat_idx]`, both in seconds of mission elapsed time, through `_maybe_add_control_tstop!`, which applies the finiteness, ordering, and span checks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `control_model` | BaseThrusterModel | n/a | yes | Positional argument `control_model`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_register_control_tstops!`; mutates `integrator` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:38-38`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:102-102`

**Downstream**

- `callees` → [[simulation.control_callbacks__maybe_add_control_tstop_bang|_maybe_add_control_tstop!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:20-20`
<!-- vulcan:connections:end -->

## Limitations
Only `start_burn_time` is bounds-checked; if `stop_burn_time` is a shorter vector the second lookup throws a `BoundsError`. The silent early return means a mis-sized model produces no tstops and no diagnostic, so the burn is integrated through rather than resolved. Nothing verifies that the stop time is later than the start time, so an inverted pair is registered as two ordinary stop points.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 16.
