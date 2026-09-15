---
id: gncx.propulsive_maneuvers_calccontroleffect_bang
label: calcControlEffect!
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: calcControlEffect!
  lines:
  - 438
  - 582
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: ControlHooks namespace dispatching the thruster control effect for
    a spacecraft index at the current simulation time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: burn_schedule
  type: PropulsiveBurnPlan
  units: mixed
  description: Validated burn window written back into the thruster model and the
    shared ODE parameter buffers.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# calcControlEffect!

## Purpose
This `calcControlEffect!` method is the thruster branch of the control contract. Each time the control callback fires for a spacecraft it decides whether a burn is currently in progress, whether one has just started or ended, and whether a new burn window should be scheduled from the commanded maneuver.

## Model & Assumptions
Burn windows come from two sources with a defined precedence: the model's own `start_burn_time` and `stop_burn_time` arrays take priority when both are finite and ordered, and the active burn plan stored in the ODE parameters is used otherwise. Once ignition begins the schedule is locked and the method returns immediately, so no replanning can occur mid-burn. After the stop time has passed, the schedule is cleared, both model times are reset to minus one, and the plan is dropped, which frees the spacecraft to plan the next campaign maneuver. Scheduling is additionally gated on flight regime: the spacecraft must be above the entry-interface altitude and pre-apoapsis, with a near-circular eccentricity tolerance that permits scheduling when the true anomaly from the orbital-element conversion is not meaningful.

## Design & Implementation
Because the control callback runs concurrently across spacecraft under `threaded_foreach_persistent(:control_callback, ...)`, the process-global maneuver-trace dictionaries are mutated only inside a `lock(_MANEUVER_TRACE_LOCK)` block. The comment in the source is explicit that Dict mutation from multiple threads is unsafe even across disjoint keys, because the internal resize and rehash are shared state. The locked block computes all four transition flags at once and returns them as a tuple, so the burn-start, burn-end, and schedule-clear trace events are emitted outside the lock. Orbital-element conversion is wrapped in a `try` that routes failures through `_control_effector_exception_fallback` rather than aborting the integration, and non-finite or hyperbolic elements cause an early return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | ControlHooks namespace dispatching the thruster control effect for a spacecraft index at the current simulation time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `burn_schedule` | PropulsiveBurnPlan | mixed | — | Validated burn window written back into the thruster model and the shared ODE parameter buffers. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:377-377`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:516-516`
- `callees` → [[core.reference_system__wrap_2pi|_wrap_2pi]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:516-516`
- `callees` → [[gnc.command_types_propulsiveburnplan|PropulsiveBurnPlan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:546-546`
- `callees` → [[gnc.propulsive_maneuvers__active_burn_plan|_active_burn_plan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:445-445`
- `callees` → [[gnc.propulsive_maneuvers__clear_burn_plan_bang|_clear_burn_plan!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:487-487`
- `callees` → [[gnc.propulsive_maneuvers__commanded_maneuver|_commanded_maneuver]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:496-496`
- `callees` → [[gnc.propulsive_maneuvers__constant_thrust_burn_duration_s|_constant_thrust_burn_duration_s]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:523-523`
- `callees` → [[gnc.propulsive_maneuvers__control_effector_exception_fallback|_control_effector_exception_fallback]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:504-504`
- `callees` → [[gnc.propulsive_maneuvers__maneuver_trace_key|_maneuver_trace_key]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:443-443`
- `callees` → [[gnc.propulsive_maneuvers__model_burn_window|_model_burn_window]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:446-446`
- `callees` → [[gnc.propulsive_maneuvers__set_burn_plan_bang|_set_burn_plan!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:558-558`
- `callees` → [[gnc.propulsive_maneuvers__trace_maneuver_event_bang|_trace_maneuver_event!]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:475-475`
- `callees` → [[gnc.propulsive_maneuvers__validated_burn_plan|_validated_burn_plan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:519-519`
- `callees` → [[parcore.reference_system_rvtoorbitalelement|rvtoorbitalelement]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:502-502`
<!-- vulcan:connections:end -->

## Limitations
Scheduling depends on osculating orbital elements computed from the instantaneous state, which are noisy near circular orbits; the eccentricity tolerance mitigates but does not remove that. The schedule lock means a maneuver cannot be aborted once started through this path. The trace dictionaries are process-global, so trace keys must be unique across concurrently running simulations in the same process.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl:438-582`.
