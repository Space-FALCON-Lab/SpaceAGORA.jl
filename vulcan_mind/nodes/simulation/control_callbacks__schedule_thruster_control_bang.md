---
id: simulation.control_callbacks__schedule_thruster_control_bang
label: _schedule_thruster_control!
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: _schedule_thruster_control!
  lines:
  - 35
  - 35
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
  description: Return value of `_schedule_thruster_control!`; mutates `integrator`
    in place. Returns `nothing`.
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

# _schedule_thruster_control!

## Purpose
Performs the full event-driven update for one thruster on one spacecraft: refresh guidance, evaluate the control effect, and register the resulting burn boundaries as solver stop points.

## Design & Implementation
Marked `@inline`, it runs three steps in fixed order. `_run_guidance_for_thruster_schedule!(integrator, sat_idx)` updates the guidance solution; `calcControlEffect!(control_model, integrator.u, integrator.p, integrator.t, sat_idx)` lets the thruster model write its thrust magnitude and the `start_burn_time`/`stop_burn_time` entries for that spacecraft; `_register_control_tstops!(integrator, control_model, sat_idx)` then feeds those two times to the integrator. The ordering matters because the tstops are read back from the model fields that `calcControlEffect!` has just written.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `control_model` | BaseThrusterModel | n/a | yes | Positional argument `control_model`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_schedule_thruster_control!`; mutates `integrator` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- [[simulation.control_callbacks_schedule_all_bang|schedule_all!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:62-62`
- [[simulation.control_callbacks_schedule_event_driven_thruster_controls_bang|schedule_event_driven_thruster_controls!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:45-45`

**Downstream**

- `callees` → [[gnc.momentum_manager_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `callees` → [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `callees` → [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `callees` → [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `callees` → [[gncx.rpo_mpc_control_model_calccontroleffect_bang|calcControlEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- `callees` → [[simulation.control_callbacks__register_control_tstops_bang|_register_control_tstops!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:38-38`
- `callees` → [[simulation.control_callbacks__run_guidance_for_thruster_schedule_bang|_run_guidance_for_thruster_schedule!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:36-36`
<!-- vulcan:connections:end -->

## Limitations
The schedule is computed from the state at the moment of the call and is not revisited unless something invokes this routine again, so a burn time that moves after scheduling leaves a stale tstop and an unresolved discontinuity. Mutation of `control_model` fields makes the routine unsafe to run concurrently for two spacecraft sharing one model unless the per-spacecraft vectors are disjoint. Failures in guidance or control propagate out of the callback.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 35.
