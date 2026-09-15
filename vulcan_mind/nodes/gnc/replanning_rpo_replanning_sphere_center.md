---
id: gnc.replanning_rpo_replanning_sphere_center
label: rpo_replanning_sphere_center
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_replanning_sphere_center
  lines:
  - 105
  - 105
inputs:
- id: sphere
  type: RPOReplanningSphere
  units: n/a
  required: true
  description: Positional argument `sphere`.
- id: t
  type: Real
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
  description: Return value of `rpo_replanning_sphere_center`. Returns `sphere.center_rtn
    + max(0.0, Float64(t) - sphere.appear_time_s) * sphere.velocit`.
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

# rpo_replanning_sphere_center

## Purpose
Computes where a drifting replanning sphere is located at simulation time `t`, applying its constant RTN velocity only after it has appeared. It is the single place that defines obstacle motion for the replanning decision and for `rpo_active_replanning_spheres`.

## Design & Implementation
Signature `rpo_replanning_sphere_center(sphere::RPOReplanningSphere, t::Real)`. It returns `sphere.center_rtn + max(0.0, Float64(t) - sphere.appear_time_s) * sphere.velocity_rtn_mps`, an `SVector{3,Float64}` in metres. The `max(0.0, ...)` clamp freezes the sphere at its declared centre for `t < appear_time_s`, so a sphere that is queried before appearing (or whose appearance time is in the future) does not move backwards along its velocity.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sphere` | RPOReplanningSphere | n/a | yes | Positional argument `sphere`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_replanning_sphere_center`. Returns `sphere.center_rtn + max(0.0, Float64(t) - sphere.appear_time_s) * sphere.velocit`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning_rpo_active_replanning_spheres|rpo_active_replanning_spheres]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:114-114`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:106-106`
<!-- vulcan:connections:end -->

## Limitations
Motion continues indefinitely after `disappear_time_s` because the function does not consult that field; callers must filter with the lifetime window first. The model is strictly linear in time, so no orbital relative motion (Clohessy-Wiltshire drift) is represented, which can misplace a real co-orbiting object over long horizons. Non-finite `t` yields a non-finite centre without an error.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 105.
