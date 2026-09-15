---
id: simulation.control_callbacks__run_guidance_for_thruster_schedule_bang
label: _run_guidance_for_thruster_schedule!
kind: function
source:
  file: src/simulation/callbacks/control_callbacks.jl
  symbol: _run_guidance_for_thruster_schedule!
  lines:
  - 28
  - 28
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
  description: Return value of `_run_guidance_for_thruster_schedule!`; mutates `integrator`
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

# _run_guidance_for_thruster_schedule!

## Purpose
Runs every configured guidance effector for spacecraft `sat_idx` immediately before the thruster schedule is computed, so the burn times the thruster model reads reflect the current guidance solution.

## Design & Implementation
Marked `@inline`. It iterates `integrator.p.args.guidance_model.guidance_effectors` under `@inbounds` and calls `calcGuidanceEffect!(guidance_model, integrator.u, integrator.p, integrator.t, sat_idx)` on each, in declaration order. Guidance effectors communicate by mutating shared state hanging off `integrator.p`, so the ordering of the vector is significant when one effector consumes another's output. The function returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `integrator` | Any | n/a | yes | Positional argument `integrator`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_run_guidance_for_thruster_schedule!`; mutates `integrator` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/control_callbacks.jl`
- [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:36-36`

**Downstream**

- `callees` → [[gnc.target_energy_bracketing_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:30-30`
- `callees` → [[gnc.thruster_guidance_functions_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:30-30`
- `callees` → [[gncy.rpo_guidance_hooks_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callers` · call · `src/simulation/callbacks/control_callbacks.jl:30-30`
<!-- vulcan:connections:end -->

## Limitations
The loop is unguarded: a single guidance effector that throws aborts the whole schedule and leaves earlier effectors' mutations applied. Order dependence is implicit, with nothing declaring which effector must run first. `@inbounds` assumes the effector vector is well formed. Because guidance is invoked here from the callback path, any effector that itself queries the atmosphere or ephemerides pays that cost inside the callback rather than in the right-hand side.

## Provenance
Mapped from `src/simulation/callbacks/control_callbacks.jl` line 28.
