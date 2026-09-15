---
id: simulation.control_callbacks__thruster_schedule_callbacks
label: _thruster_schedule_callbacks
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: _thruster_schedule_callbacks
  lines:
  - 51
  - 51
inputs:
- id: control_model
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `control_model`.
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
  type: Union{Nothing, Tuple}
  units: n/a
  description: Return value of `_thruster_schedule_callbacks`. Returns `nothing` or
    `(init_callback,)`.
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

# _thruster_schedule_callbacks

## Purpose
Builds the callback tuple that installs an event-driven thruster schedule for all `num_sats` spacecraft once, at solver initialization, instead of polling the thruster at a fixed rate.

## Design & Implementation
It first validates shape: `length(control_model.thrust)` must equal `num_sats`, otherwise it throws an `ArgumentError` naming both counts and directing the user to a single shared model carrying per-spacecraft vectors. It then defines the closure `schedule_all!`, which loops `sat_idx in 1:num_sats` calling `_schedule_thruster_control!`, and wraps it in a `DiscreteCallback` whose condition is the constant `(u, t, integrator) -> false` so it never fires during the run, while its `initialize` hook calls `schedule_all!(integrator)`. The result is returned as a one-element tuple so the caller can `append!` it uniformly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `control_model` | BaseThrusterModel | n/a | yes | Positional argument `control_model`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, Tuple} | n/a | — | Return value of `_thruster_schedule_callbacks`. Returns `nothing` or `(init_callback,)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:93-93`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the condition is permanently false, the schedule is computed exactly once from the initial state; any guidance update later in the run does not refresh the burn times through this path and must call `schedule_event_driven_thruster_controls!` explicitly. The length check covers `thrust` only, so mismatched `start_burn_time` or `stop_burn_time` vectors are caught later as bounds errors. The whole spacecraft loop runs serially inside the initialization hook.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 51.
