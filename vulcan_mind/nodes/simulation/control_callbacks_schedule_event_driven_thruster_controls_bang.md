---
id: simulation.control_callbacks_schedule_event_driven_thruster_controls_bang
label: schedule_event_driven_thruster_controls!
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: schedule_event_driven_thruster_controls!
  lines:
  - 42
  - 42
inputs:
- id: integrator
  type: Any
  units: n/a
  required: true
  description: Positional argument `integrator`.
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
  description: Return value of `schedule_event_driven_thruster_controls!`; mutates
    `integrator` in place. Returns `nothing`.
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

# schedule_event_driven_thruster_controls!

## Purpose
Rebuilds the event-driven thruster schedule for a single spacecraft across every thruster effector in the control model. Used when an external event, such as a replan, invalidates the burn times computed at initialization.

## Design & Implementation
It walks `integrator.p.args.control_model.control_effectors` under `@inbounds`, tests each entry with `control_model isa BaseThrusterModel`, and calls `_schedule_thruster_control!(integrator, control_model, sat_idx)` for the matches. Non-thruster effectors are skipped because they are driven by their own `PeriodicCallback` at the configured control rate. The function targets one `sat_idx` at a time, which keeps it usable from inside a per-spacecraft threaded loop as well as from single-spacecraft event handlers.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `schedule_event_driven_thruster_controls!`; mutates `integrator` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- [[simulation.event_callbacks_affect_upcrossing_bang|affect_upcrossing!]] · `callees` → `callers` · call · `src/simulation/callbacks/event_callbacks.jl:169-169`

**Downstream**

- `callees` → [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:45-45`
<!-- vulcan:connections:end -->

## Limitations
Re-running this adds new tstops without removing previously registered ones, so a repeatedly rescheduled burn leaves obsolete stop points that cost extra solver steps. The `isa` filter means a thruster wrapped in a composite or decorator effector is not recognised and never gets scheduled. As with the single-thruster path, it mutates per-spacecraft entries of shared effector objects.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 42.
