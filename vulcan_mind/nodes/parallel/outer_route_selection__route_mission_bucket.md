---
id: parallel.outer_route_selection__route_mission_bucket
label: _route_mission_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_mission_bucket
  lines:
  - 40
  - 40
inputs:
- id: mission_time_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_time_s`.
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
  type: String
  units: n/a
  description: Return value of `_route_mission_bucket`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _route_mission_bucket

## Purpose
Classifies mission duration in seconds into "short", "medium", or "long" for the `mission=` signature field, so routing feedback is shared among simulations of comparable wall-clock scale.

## Design & Implementation
`@inline _route_mission_bucket(mission_time_s::Float64)::String` returns "short" for `<= 1800.0` s (30 min), "medium" for `<= 7200.0` s (2 h), and "long" above that. Used by every signature builder in the file. Note that `_feature_is_lightweight` and `_feature_heavy_for_process` use the separately tunable `t.outer_light_mission_threshold_s` rather than these constants.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mission_time_s` | Float64 | n/a | yes | Positional argument `mission_time_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_mission_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__compat_outer_route_signature|_compat_outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:139-139`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:188-188`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:160-160`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Thresholds are hard-coded and unrelated to the tuning thresholds, so a mission can be "long" for signature purposes yet lightweight for routing, or vice versa. `NaN` compares false on both branches and lands in "long"; negative durations become "short".

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 40.
