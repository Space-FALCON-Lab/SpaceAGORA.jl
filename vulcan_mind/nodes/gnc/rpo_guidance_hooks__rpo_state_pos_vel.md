---
id: gnc.rpo_guidance_hooks__rpo_state_pos_vel
label: _rpo_state_pos_vel
kind: function
source:
  file: src/gnc/guidance/rpo/rpo_guidance_hooks.jl
  symbol: _rpo_state_pos_vel
  lines:
  - 2
  - 2
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  type: SVector
  units: n/a
  description: Return value of `_rpo_state_pos_vel`. Returns `SVector{3, Float64}(u.sc[idx].pos),
    SVector{3, Float64}(u.sc[idx].vel)`.
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

# _rpo_state_pos_vel

## Purpose
Pulls the inertial position and velocity of one spacecraft out of the packed simulation state `u`, converting them to stack-allocated static vectors for the RPO relative-motion computations.

## Design & Implementation
Takes the state container `u` and an integer `idx`, and returns the tuple `(SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel))`. The explicit `SVector{3, Float64}` construction both fixes the element type, so downstream arithmetic stays type-stable even when the integrator is running a dual-number or extended-precision state, and avoids aliasing the mutable state arrays. Both `build_rpo_plan` and `maybe_update_rpo_replanning!` call it twice per invocation, once for `model.chaser_idx` and once for `model.target_idx`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `_rpo_state_pos_vel`. Returns `SVector{3, Float64}(u.sc[idx].pos), SVector{3, Float64}(u.sc[idx].vel)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_build_rpo_plan|build_rpo_plan]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:9-9`
- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:89-89`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The index is not range-checked against the number of spacecraft in `u.sc`, so a misconfigured `chaser_idx` or `target_idx` surfaces as a raw `BoundsError` with no guidance-level context. Because the result is a copy, later mutation of the state is not reflected; callers must re-read every step. The conversion to `Float64` silently drops derivative information if `u` carries dual numbers for sensitivity analysis.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_guidance_hooks.jl` line 2.
