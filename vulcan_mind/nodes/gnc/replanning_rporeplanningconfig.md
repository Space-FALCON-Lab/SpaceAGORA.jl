---
id: gnc.replanning_rporeplanningconfig
label: RPOReplanningConfig
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: RPOReplanningConfig
  lines:
  - 29
  - 29
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: true
  description: Field `enabled`.
- id: spheres
  type: Vector{RPOReplanningSphere}
  units: n/a
  required: true
  description: Field `spheres`.
- id: desired_goal_rtn
  type: Union{Nothing, SVector{3, Float64}}
  units: n/a
  required: true
  description: Field `desired_goal_rtn`.
- id: goal_change_tolerance_m
  type: Float64
  units: n/a
  required: true
  description: Field `goal_change_tolerance_m`.
- id: safe_distance_m
  type: Float64
  units: n/a
  required: true
  description: Field `safe_distance_m`.
- id: retime_clearance_m
  type: Float64
  units: n/a
  required: true
  description: Field `retime_clearance_m`.
- id: min_replan_interval_s
  type: Float64
  units: n/a
  required: true
  description: Field `min_replan_interval_s`.
- id: hysteresis_samples
  type: Int
  units: n/a
  required: true
  description: Field `hysteresis_samples`.
- id: sphere_surface_samples
  type: Int
  units: n/a
  required: true
  description: Field `sphere_surface_samples`.
- id: remaining_sample_ds_m
  type: Float64
  units: n/a
  required: true
  description: Field `remaining_sample_ds_m`.
- id: tracking_error_retime_m
  type: Float64
  units: n/a
  required: true
  description: Field `tracking_error_retime_m`.
- id: tracking_error_replan_m
  type: Float64
  units: n/a
  required: true
  description: Field `tracking_error_replan_m`.
- id: rng_seed
  type: Int
  units: n/a
  required: true
  description: Field `rng_seed`.
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
  type: RPOReplanningConfig
  units: n/a
  description: Constructed `RPOReplanningConfig`.
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

# RPOReplanningConfig

## Purpose
Immutable policy object that tells the RPO guidance loop when to keep, retime or replan the active HYPR path in response to dynamic obstacles, a changed goal, or accumulated tracking error. It bundles the obstacle list with every threshold the decision logic in `rpo_replanning_decision` consults.

## Design & Implementation
A plain `struct RPOReplanningConfig` with fields `enabled::Bool`, `spheres::Vector{RPOReplanningSphere}`, `desired_goal_rtn::Union{Nothing,SVector{3,Float64}}`, `goal_change_tolerance_m`, `safe_distance_m`, `retime_clearance_m`, `min_replan_interval_s`, `hysteresis_samples::Int`, `sphere_surface_samples::Int`, `remaining_sample_ds_m`, `tracking_error_retime_m`, `tracking_error_replan_m` and `rng_seed::Int`. The keyword constructor (defaults: `enabled=true`, `goal_change_tolerance_m=1e-6`, `safe_distance_m=0.0`, `retime_clearance_m=0.25`, `min_replan_interval_s=0.0`, `hysteresis_samples=1`, `sphere_surface_samples=96`, `remaining_sample_ds_m=0.10`, both tracking thresholds `Inf`, `rng_seed=740`) normalises each sphere-like input through `_rpo_replanning_sphere`, clamps the tracking thresholds at zero, and throws `ArgumentError` for negative tolerances, `retime_clearance_m < safe_distance_m`, `hysteresis_samples < 1`, `sphere_surface_samples < 12`, non-positive `remaining_sample_ds_m`, or NaN tracking thresholds.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | yes | Field `enabled`. |
| in | `spheres` | Vector{RPOReplanningSphere} | n/a | yes | Field `spheres`. |
| in | `desired_goal_rtn` | Union{Nothing, SVector{3, Float64}} | n/a | yes | Field `desired_goal_rtn`. |
| in | `goal_change_tolerance_m` | Float64 | n/a | yes | Field `goal_change_tolerance_m`. |
| in | `safe_distance_m` | Float64 | n/a | yes | Field `safe_distance_m`. |
| in | `retime_clearance_m` | Float64 | n/a | yes | Field `retime_clearance_m`. |
| in | `min_replan_interval_s` | Float64 | n/a | yes | Field `min_replan_interval_s`. |
| in | `hysteresis_samples` | Int | n/a | yes | Field `hysteresis_samples`. |
| in | `sphere_surface_samples` | Int | n/a | yes | Field `sphere_surface_samples`. |
| in | `remaining_sample_ds_m` | Float64 | n/a | yes | Field `remaining_sample_ds_m`. |
| in | `tracking_error_retime_m` | Float64 | n/a | yes | Field `tracking_error_retime_m`. |
| in | `tracking_error_replan_m` | Float64 | n/a | yes | Field `tracking_error_replan_m`. |
| in | `rng_seed` | Int | n/a | yes | Field `rng_seed`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOReplanningConfig | n/a | — | Constructed `RPOReplanningConfig`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning__rpo_replanning_sphere|_rpo_replanning_sphere]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:68-68`
- [[gnc.rpo_guidance_hooks__rpo_replanning_config|_rpo_replanning_config]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:63-63`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`min_replan_interval_s`, `hysteresis_samples` and `rng_seed` are validated and stored but not consumed anywhere in this file; enforcement, if any, must happen in the caller. Negative tracking thresholds are silently clamped to zero rather than rejected, unlike the other fields. The default constructor without keywords is still callable positionally and skips all checks. `spheres` is a mutable `Vector`, so the config is not deeply immutable.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 29.
