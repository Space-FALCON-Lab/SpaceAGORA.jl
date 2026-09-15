---
id: gncy.replanning_rpo_replanning_decision
label: rpo_replanning_decision
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_replanning_decision
  lines:
  - 211
  - 246
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: plan_state
  type: Tuple
  units: n/a
  required: true
  description: Active RPO plan, current RTN position, station geometry, replanning
    configuration, and the current mission time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: decision
  type: NamedTuple
  units: n/a
  description: Action symbol, reason symbol, minimum clearance, active spheres, and
    the geometry the action should be executed against.
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

# rpo_replanning_decision

## Purpose
`rpo_replanning_decision` is the in-flight supervisor for an active RPO plan. On each guidance update it decides whether the current reference remains acceptable, needs only a new time parameterisation, or must be replanned geometrically, and it returns the augmented geometry the chosen action should use.

## Model & Assumptions
The decision uses a strict priority order rather than a scored blend. A goal change outranks everything, because a plan aimed at the wrong terminal state cannot be salvaged by retiming. Dynamic obstacles outrank tracking error, because a keep-out intrusion is a safety condition while tracking error is a performance condition. Within each branch, replanning outranks retiming, so the more expensive but more capable action is chosen whenever the cheaper one would be insufficient. A disabled configuration or an invalid plan short-circuits to no action.

## Design & Implementation
Active dynamic obstacles are gathered by `rpo_active_replanning_spheres` at the current time and folded into the static station geometry by `rpo_geometry_with_replanning_spheres` using `sphere_surface_samples` points per sphere. When a desired goal is configured, the last column of the reference position matrix is compared against it and a mismatch beyond `goal_change_tolerance_m` triggers an immediate replan. With no active spheres the routine falls back to `rpo_reference_tracking_error`, a brute-force nearest-reference-point search, and compares it against `tracking_error_replan_m` and then `tracking_error_retime_m`. With active spheres it extracts the remaining reference path from the current position at `remaining_sample_ds_m` spacing and evaluates clearance statistics against the augmented geometry, replanning below `safe_distance_m` and retiming below `retime_clearance_m`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `plan_state` | Tuple | n/a | yes | Active RPO plan, current RTN position, station geometry, replanning configuration, and the current mission time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `decision` | NamedTuple | n/a | — | Action symbol, reason symbol, minimum clearance, active spheres, and the geometry the action should be executed against. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rpo_guidance_hooks_maybe_update_rpo_replanning_bang|maybe_update_rpo_replanning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:93-93`

**Downstream**

- `callees` → [[gnc.clearance_rpo_path_clearance_stats|rpo_path_clearance_stats]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:239-239`
- `callees` → [[gnc.replanning_rpo_active_replanning_spheres|rpo_active_replanning_spheres]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:215-215`
- `callees` → [[gnc.replanning_rpo_geometry_with_replanning_spheres|rpo_geometry_with_replanning_spheres]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:216-216`
- `callees` → [[gnc.replanning_rpo_reference_tracking_error|rpo_reference_tracking_error]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:229-229`
- `callees` → [[gnc.replanning_rpo_remaining_reference_path|rpo_remaining_reference_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:238-238`
<!-- vulcan:connections:end -->

## Limitations
Tracking error is measured as the distance to the nearest reference point over the entire plan, so a vehicle that is close to a much earlier or later part of the reference reports a small error even when it is badly off schedule. The goal comparison inspects only the final reference column and so does not detect a goal that moved within tolerance repeatedly. Tracking-error thresholds are ignored whenever any sphere is active, so a large tracking error during a dynamic-obstacle episode does not by itself trigger action. Clearance is evaluated only on the remaining path, so an obstacle that will appear later is invisible until its appear time.

## Provenance
Mapped from replanning.jl lines 197-246; include site observed at guidance_hooks.jl line 77.
