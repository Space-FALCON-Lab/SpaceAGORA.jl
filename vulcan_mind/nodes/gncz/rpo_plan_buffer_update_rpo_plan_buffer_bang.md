---
id: gncz.rpo_plan_buffer_update_rpo_plan_buffer_bang
label: update_rpo_plan_buffer!
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_plan_buffer.jl
  symbol: update_rpo_plan_buffer!
  lines:
  - 20
  - 25
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: GNC guidance namespace holding the RPO plan record types and the buffer
    update entry point.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: buffer
  type: RPOPlanBuffer
  units: n/a
  description: The same buffer object after the new plan, its validity flag, and the
    update timestamp have been written in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# update_rpo_plan_buffer!

## Purpose
`update_rpo_plan_buffer!` installs a freshly computed rendezvous and proximity operations plan into the mutable buffer that guidance reads during propagation. The buffer is the single place where a planner result becomes visible to the rest of the simulation, so the update has to leave the buffer in a state that a consumer can interpret without knowing whether planning succeeded.

## Model & Assumptions
The file declares two keyword-constructed mutable records. `RPOPlan` carries a validity flag, a reference time vector, reference position and velocity matrices in the radial-transverse-normal frame, the raw geometric path, a scalar cost, and a named tuple of planner diagnostics. Every field defaults to an empty or infinite value, so a default-constructed plan represents the absence of a usable trajectory rather than a zero-length trajectory. `RPOPlanBuffer` wraps one plan together with its own validity flag and the simulation time at which it was installed, defaulting that timestamp to a not-a-number sentinel.

## Design & Implementation
The update assigns the plan, copies the plan validity onto the buffer so a consumer can gate on one field, converts the supplied time to a double precision value, and returns the buffer to allow chaining. Buffer validity is therefore always derived from the plan rather than being set independently, which prevents a stale valid flag from surviving an installation of an invalid plan.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | GNC guidance namespace holding the RPO plan record types and the buffer update entry point. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `buffer` | RPOPlanBuffer | n/a | — | The same buffer object after the new plan, its validity flag, and the update timestamp have been written in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:116-116`
- [[gncy.rpo_guidance_hooks_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:155-155`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/rpo_plan_buffer.jl:23-23`
<!-- vulcan:connections:end -->

## Limitations
The buffer stores only the current plan, so the docstring promise of retaining a previous plan for diagnostics is not realised by the assignment. Nothing checks that the reference matrices are dimensionally consistent, that the time vector matches the number of reference columns, or that the installation time moves forward. Because the plan is stored by reference, later mutation of the caller's plan object is visible through the buffer.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_plan_buffer.jl:1-25`.
