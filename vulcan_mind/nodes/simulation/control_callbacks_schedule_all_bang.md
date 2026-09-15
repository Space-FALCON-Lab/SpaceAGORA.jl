---
id: simulation.control_callbacks_schedule_all_bang
label: schedule_all!
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: schedule_all!
  lines:
  - 60
  - 60
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
  type: Nothing
  units: n/a
  description: Return value of `schedule_all!`; mutates `integrator` in place. Returns
    `nothing`.
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

# schedule_all!

## Purpose
Initialization body of the thruster schedule callback. It walks every spacecraft from 1 to `num_sats` and lays down that spacecraft's burn schedule and solver stop points in one pass.

## Design & Implementation
Defined as a closure inside `_thruster_schedule_callbacks`, capturing `control_model` and `num_sats`, so it can serve both as the callback's `affect!` and as its `initialize` hook without re-deriving state. The body is a single `@inbounds for sat_idx in 1:num_sats` loop over `_schedule_thruster_control!(integrator, control_model, sat_idx)`, which runs guidance, evaluates the control effect, and registers the start and stop burn tstops for that spacecraft. It returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `schedule_all!`; mutates `integrator` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`

**Downstream**

- `callees` → [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:62-62`
<!-- vulcan:connections:end -->

## Limitations
The loop is strictly serial, so initialization cost grows linearly with the constellation size even when the surrounding simulation is threaded. There is no error isolation between spacecraft: a guidance or control failure on one index aborts scheduling for every later index, leaving a partially populated tstop set. Since the enclosing callback condition is always false, this closure normally executes only once at `t0`.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 60.
