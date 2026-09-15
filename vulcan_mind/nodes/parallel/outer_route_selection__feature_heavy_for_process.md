---
id: parallel.outer_route_selection__feature_heavy_for_process
label: _feature_heavy_for_process
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _feature_heavy_for_process
  lines:
  - 274
  - 274
inputs:
- id: f
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `f`.
- id: t
  type: OuterRouteTuning
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
  type: Bool
  units: n/a
  description: Return value of `_feature_heavy_for_process`.
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

# _feature_heavy_for_process

## Purpose
Decides whether a non-Monte-Carlo workload is expensive enough to justify spawning worker processes, which gates whether `:process` is added to the candidate list in `outer_route_candidates`.

## Design & Implementation
`@inline _feature_heavy_for_process(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool` returns true immediately for native GRAM point density (lock-limited, comment in source), false if `_feature_is_lightweight(f, t)`, and otherwise `f.has_nbody || f.harmonics_degree >= 20 || f.mission_time_s > t.outer_light_mission_threshold_s`. The mission-time clause means any workload longer than the lightweight threshold that also fails lightweight for another reason is process-eligible.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `t` | OuterRouteTuning | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_feature_heavy_for_process`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection_outer_route_candidates|outer_route_candidates]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:398-398`

**Downstream**

- `callees` → [[parallel.outer_route_selection__feature_is_lightweight|_feature_is_lightweight]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:279-279`
- `callees` → [[parallel.outer_route_selection__is_native_gram_point_density|_is_native_gram_point_density]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:275-275`
<!-- vulcan:connections:end -->

## Limitations
The harmonics threshold of 20 is a literal, separate from any tuning field. `machine_class` is not consulted here, so a small machine can still be offered `:process` as a candidate; only `default_outer_route` weighs machine class. SRP, thermal, and effector cost are not considered.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 274.
