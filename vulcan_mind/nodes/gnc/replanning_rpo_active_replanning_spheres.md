---
id: gnc.replanning_rpo_active_replanning_spheres
label: rpo_active_replanning_spheres
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_active_replanning_spheres
  lines:
  - 110
  - 110
inputs:
- id: config
  type: RPOReplanningConfig
  units: n/a
  required: true
  description: Positional argument `config`.
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
  description: Return value of `rpo_active_replanning_spheres`. Returns `RPOReplanningSphere[`.
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

# rpo_active_replanning_spheres

## Purpose
Produces the list of obstacles that exist at time `t`, each re-centred at its drifted position, so that downstream geometry augmentation and clearance checks can treat them as static spheres for that instant.

## Design & Implementation
Signature `rpo_active_replanning_spheres(config::RPOReplanningConfig, t::Real)`, returning a `Vector{RPOReplanningSphere}`. It converts `t` to `Float64`, iterates `config.spheres`, keeps those with `appear_time_s <= t <= disappear_time_s`, and for each constructs a new `RPOReplanningSphere` via the validating constructor with `center_rtn = rpo_replanning_sphere_center(sphere, t)` and the original `radius_m`, `appear_time_s`, `disappear_time_s`, `velocity_rtn_mps` and `label`. The original config is not mutated.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | RPOReplanningConfig | n/a | yes | Positional argument `config`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_active_replanning_spheres`. Returns `RPOReplanningSphere[`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.replanning_rpo_replanning_decision|rpo_replanning_decision]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:215-215`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:111-111`
- `callees` → [[gnc.replanning_rpo_replanning_sphere_center|rpo_replanning_sphere_center]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:114-114`
- `callees` → [[gnc.replanning_rporeplanningsphere|RPOReplanningSphere]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:113-113`
<!-- vulcan:connections:end -->

## Limitations
The returned spheres retain their `velocity_rtn_mps` and `appear_time_s`, so calling `rpo_replanning_sphere_center` on them again would apply the drift a second time; they are only correct as instantaneous snapshots. Each call allocates a new vector and re-runs constructor validation for every active sphere, which is wasteful inside a high-rate guidance loop. Inclusive bounds mean a sphere is active at exactly its disappearance instant.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 110.
