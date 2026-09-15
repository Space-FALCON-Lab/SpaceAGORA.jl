---
id: parallel.outer_route_selection__feature_is_lightweight
label: _feature_is_lightweight
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _feature_is_lightweight
  lines:
  - 258
  - 258
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
  description: Return value of `_feature_is_lightweight`.
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

# _feature_is_lightweight

## Purpose
Predicate that identifies workloads cheap enough that any outer parallel route would cost more in overhead than it saves, so `default_outer_route` returns `:none` and `_feature_heavy_for_process` returns `false` for them.

## Design & Implementation
`@inline _feature_is_lightweight(f::OuterRouteFeatures, t::OuterRouteTuning)::Bool` is true only when `f.n_sats <= t.outer_light_sat_threshold`, `f.n_links <= t.outer_light_link_threshold`, `f.mission_time_s <= t.outer_light_mission_threshold_s`, and the workload has no n-body gravity, no control, no orientation dynamics, and `harmonics_degree == 0`. All thresholds come from the tuning struct, making this the one routing rule that is fully configurable.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `t` | OuterRouteTuning | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_feature_is_lightweight`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.default_outer_route|default_outer_route]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:337-337`
- [[parallel.outer_route_selection__feature_heavy_for_process|_feature_heavy_for_process]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:279-279`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Atmosphere model, SRP, thermal, and effector counts are ignored, so a single-satellite native-GRAM run with heavy effectors can still be labelled lightweight here (the native-GRAM case is caught earlier by `_is_native_gram_point_density`, but expensive non-GRAM density models are not). Monte Carlo sample count is also not considered.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 258.
