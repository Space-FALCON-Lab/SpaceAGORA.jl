---
id: gnc.rpo_guidance_hooks_build_rpo_plan
label: build_rpo_plan
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_guidance_hooks.jl
  symbol: build_rpo_plan
  lines:
  - 7
  - 7
inputs:
- id: model
  type: RPOGuidanceModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Any
  units: n/a
  description: Return value of `build_rpo_plan`. Returns `build_rpo_plan_from_start(model,
    x_rel[1:3], model.geometry, t)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# build_rpo_plan

## Purpose
Builds a complete RPO trajectory plan from scratch using the guidance model's stored goal and the live relative state of chaser and target at simulation time `t`.

## Design & Implementation
First guards the reference geometry: `model.geometry === nothing` throws `ArgumentError("RPOGuidanceModel requires reference geometry.")`. It then reads both spacecraft states through `_rpo_state_pos_vel`, converts them into the target-centred RTN frame with `inertial_to_rtn_relative_state(r_chaser, v_chaser, r_target, v_target)`, and passes only the position part, `x_rel[1:3]`, along with `model.geometry` and `t` to `build_rpo_plan_from_start`. The relative velocity components of `x_rel` are computed but deliberately discarded, because the planner searches in position space and reconstructs a velocity reference afterwards.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | RPOGuidanceModel | n/a | yes | Positional argument `model`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `build_rpo_plan`. Returns `build_rpo_plan_from_start(model, x_rel[1:3], model.geometry, t)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.rpo_guidance_hooks_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:154-154`

**Downstream**

- `callees` → [[core.reference_system_inertial_to_rtn_relative_state|inertial_to_rtn_relative_state]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:11-11`
- `callees` → [[gnc.rpo_guidance_hooks__rpo_state_pos_vel|_rpo_state_pos_vel]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:9-9`
- `callees` → [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:12-12`
<!-- vulcan:connections:end -->

## Limitations
Discarding the relative velocity means the resulting path does not start tangent to the chaser's current motion, so the reference can demand a discontinuous velocity change at the handover instant. The geometry snapshot is taken by reference from the model rather than copied, so a dynamic obstacle set that mutates during planning is seen inconsistently. No safe-distance override is offered on this entry point, and there is no timeout, so a difficult planning problem blocks the simulation step for as long as the swarm optimiser runs.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_guidance_hooks.jl` line 7.
