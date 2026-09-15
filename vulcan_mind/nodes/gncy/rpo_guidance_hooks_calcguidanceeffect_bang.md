---
id: gncy.rpo_guidance_hooks_calcguidanceeffect_bang
label: calcGuidanceEffect!
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_guidance_hooks.jl
  symbol: calcGuidanceEffect!
  lines:
  - 151
  - 161
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: sim_state
  type: Tuple
  units: n/a
  required: true
  description: Guidance model, simulation state vector, ODE parameters, current time,
    and the index of the spacecraft being updated.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: plan_buffer_update
  type: Nothing
  units: n/a
  description: In-place update of the guidance model plan buffer and replanning counters;
    the call returns nothing.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# calcGuidanceEffect!

## Purpose
`calcGuidanceEffect!` is the per-step guidance entry point the simulation loop calls for an RPO guidance model. It decides between building a fresh plan and asking the replanning supervisor whether the existing plan still holds, and it writes the outcome into the model plan buffer.

## Model & Assumptions
Only the chaser is planned for. The function returns immediately unless the spacecraft index matches `model.chaser_idx`, so the target vehicle is propagated without guidance even though both appear in the same state vector. Planning is triggered either by an invalid buffer, which is the cold-start case, or by an explicit `force_replan` flag, which lets mission logic demand a fresh plan without inspecting buffer internals. Once a valid plan exists the model switches permanently into supervised mode and all further changes flow through `maybe_update_rpo_replanning!`.

## Design & Implementation
The cold path calls `build_rpo_plan`, which extracts chaser and target position and velocity through `_rpo_state_pos_vel`, converts to a relative RTN state with `inertial_to_rtn_relative_state`, and forwards the position part to `build_rpo_plan_from_start`. That function resolves the safe distance from the model or the configuration, runs `rpo_pso_plan_path`, converts the geometric path to a timed reference with `rpo_reference_from_path`, and packs an `RPOPlan` whose diagnostics carry the cost components, adaptive record, refinement flag, early-stop status, iteration timeout status and events, warm-start diagnostics, and the planning time. The result is written by `update_rpo_plan_buffer!` and `force_replan` is cleared. Replanning failures are caught in `maybe_update_rpo_replanning!`, which restores the original goal, increments `replan_failure_count`, records a `:replan_failed` event, and warns while keeping the previous plan.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `sim_state` | Tuple | n/a | yes | Guidance model, simulation state vector, ODE parameters, current time, and the index of the spacecraft being updated. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `plan_buffer_update` | Nothing | n/a | — | In-place update of the guidance model plan buffer and replanning counters; the call returns nothing. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.control_callbacks__run_guidance_for_thruster_schedule_bang|_run_guidance_for_thruster_schedule!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:30-30`
- [[simulation.navigation_guidance_callbacks_get_guidance_callbacks|get_guidance_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/navigation_guidance_callbacks.jl:27-27`

**Downstream**

- `callees` → [[gnc.rpo_guidance_hooks_build_rpo_plan|build_rpo_plan]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:154-154`
- `callees` → [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:158-158`
- `callees` → [[gncz.rpo_plan_buffer_update_rpo_plan_buffer_bang|update_rpo_plan_buffer!]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:155-155`
<!-- vulcan:connections:end -->

## Limitations
Planning happens inline on the guidance step, so a cold start or a forced replan stalls the integrator for the full duration of a PSO solve. Missing geometry is reported only as an argument error thrown from deep inside `build_rpo_plan`. The chaser index guard means a mis-configured index silently produces no guidance at all rather than an error. Replanning failures preserve the previous plan even when the supervisor judged it unsafe, so a persistently failing planner leaves the vehicle following a reference that already triggered a replan request.

## Provenance
Mapped from rpo_guidance_hooks.jl lines 1-161; include site observed at guidance_hooks.jl line 78.
